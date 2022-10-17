#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

uint8_t field_mul(uint8_t a, uint8_t b) {
    uint8_t log[] = {0, 15, 4, 1, 8, 2, 5, 10, 12, 11, 6, 13, 9, 7, 14, 3};
    uint8_t exp[] = {1, 3, 5, 15, 2, 6, 10, 13, 4, 12, 7, 9, 8, 11, 14};
    if (a == 0 || b == 0)
        return 0;
    return exp[(log[a] + log[b]) % 15];
}

uint8_t field_inv(uint8_t a) {
    uint8_t t[] = {0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8};
    return t[a];
}

void erasure_code_build_table(
    uint8_t *table,
    int input_count,
    int output_count,
    const int *input_x,
    const int *output_x
) {
    assert(input_count + output_count <= 16);
    for (int i = 0; i < input_count; i++) {
        assert(0 <= input_x[i] && input_x[i] < 16);
        for (int j = 0; j < output_count; j++) {
            assert(0 <= output_x[j] && output_x[j] < 16);
            uint8_t accum = 1;
            for (int p = 0; p < input_count; p++) {
                if (p == i)
                    continue;
                accum = field_mul(accum, field_mul(
                    output_x[j] ^ input_x[p],
                    field_inv(input_x[i] ^ input_x[p])
                ));
            }
            for (int bit = 0; bit < 4; bit++) {
                uint8_t output_pattern = field_mul(accum, 1 << bit);
                table[i * output_count * 4 + j * 4 + bit] = output_pattern;
            }
        }
    }
}

void erasure_code_process_old(
    uint8_t *table,
    int input_count,
    int output_count,
    const uint8_t **input_shard_data,
    uint8_t **output_shard_data,
    size_t shard_length
) {
    assert(shard_length % 32 == 0);
    const uint64_t **input64 = (const uint64_t **)input_shard_data;
    uint64_t **output64 = (uint64_t **)output_shard_data;
    size_t segments = shard_length / 32;
    for (int i = 0; i < input_count; i++) {
        for (int j = 0; j < output_count; j++) {
            for (int segment = 0; segment < segments; segment++) {
                const uint64_t *input_ptr = input64[i] + segment * 4;
                uint64_t *output_ptr = output64[j] + segment * 4;
                for (int bit_in = 0; bit_in < 4; bit_in++) {
                    uint8_t mask = table[i * output_count * 4 + j * 4 + bit_in];
                    if (mask & 1) output_ptr[0] ^= input_ptr[bit_in];
                    if (mask & 2) output_ptr[1] ^= input_ptr[bit_in];
                    if (mask & 4) output_ptr[2] ^= input_ptr[bit_in];
                    if (mask & 8) output_ptr[3] ^= input_ptr[bit_in];
                }
            }
        }
    }
}

void erasure_code_process(
    uint8_t *table,
    int input_count,
    int output_count,
    const uint8_t **input_shard_data,
    uint8_t **output_shard_data,
    size_t shard_length
) {
    assert(shard_length % 24 == 0);
    const uint64_t **data = (const uint64_t **)input_shard_data;
    volatile uint64_t **out = (uint64_t **)output_shard_data;
    size_t segments = shard_length / 24;
    for (int i = 0; i < segments; i++) {
        out[0][3*i + 0] = data[0][3*i + 0];// ^ data[3][3*i + 0];
        out[0][3*i + 1] = data[0][3*i + 1];// ^ data[3][3*i + 1];
        out[0][3*i + 2] = data[0][3*i + 2];// ^ data[3][3*i + 2];

        out[1][3*i + 0] = data[1][3*i + 0];// ^ data[4][3*i + 0];
        out[1][3*i + 1] = data[1][3*i + 1];// ^ data[4][3*i + 1];
        out[1][3*i + 2] = data[1][3*i + 2];// ^ data[4][3*i + 2];

        out[2][3*i + 0] = data[2][3*i + 0];// ^ data[0][3*i + 0];
        out[2][3*i + 1] = data[2][3*i + 1];// ^ data[0][3*i + 1];
        out[2][3*i + 2] = data[2][3*i + 2];// ^ data[0][3*i + 2];

        // out[0][3*i + 0] = data[0][3*i + 0] ^ data[0][3*i + 1] ^ data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[1][3*i + 2] ^ data[2][3*i + 1] ^ data[2][3*i + 2] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 0];
        // out[0][3*i + 1] = data[0][3*i + 0] ^ data[1][3*i + 0] ^ data[2][3*i + 0] ^ data[2][3*i + 1] ^ data[3][3*i + 0] ^ data[3][3*i + 1] ^ data[4][3*i + 1];
        // out[0][3*i + 2] = data[0][3*i + 0] ^ data[0][3*i + 1] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[2][3*i + 0] ^ data[2][3*i + 1] ^ data[2][3*i + 2] ^ data[3][3*i + 0] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 2];

        // out[1][3*i + 0] = data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[1][3*i + 2] ^ data[2][3*i + 2] ^ data[3][3*i + 0] ^ data[3][3*i + 2] ^ data[4][3*i + 0];
        // out[1][3*i + 1] = data[0][3*i + 0] ^ data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[1][3*i + 2] ^ data[2][3*i + 0] ^ data[2][3*i + 2] ^ data[3][3*i + 0] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 1];
        // out[1][3*i + 2] = data[0][3*i + 1] ^ data[1][3*i + 1] ^ data[1][3*i + 2] ^ data[2][3*i + 1] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 2];

        // out[2][3*i + 0] = data[0][3*i + 1] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[2][3*i + 0] ^ data[2][3*i + 1] ^ data[3][3*i + 1] ^ data[4][3*i + 0];
        // out[2][3*i + 1] = data[0][3*i + 1] ^ data[0][3*i + 2] ^ data[1][3*i + 2] ^ data[2][3*i + 2] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 1];
        // out[2][3*i + 2] = data[0][3*i + 0] ^ data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[2][3*i + 0] ^ data[3][3*i + 0] ^ data[3][3*i + 2] ^ data[4][3*i + 2];
    }
}


#define BYTES (2 * 1024 * 1024 + 16)

int main() {
    // Configure to encode 5 data shards into 3 parity shards.
    uint8_t table[5 * 3 * 4];
    int input_x[5] = {0, 1, 2, 3, 4};
    int output_x[3] = {5, 6, 7};
    erasure_code_build_table(table, 5, 3, input_x, output_x);

    uint8_t *shards[8];
    for (int i = 0; i < 8; i++) {
        shards[i] = malloc(BYTES);
        for (int j = 0; j < BYTES; j++) {
            shards[i][j] = i * BYTES + j;
        }
    }

    for (int j = 0; j < 1000; j++) {
        erasure_code_process(table, 5, 3, (const uint8_t**)shards, shards + 5, BYTES);
    }
}
