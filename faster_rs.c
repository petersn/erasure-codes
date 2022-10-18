#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <x86intrin.h>

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


__attribute__((noinline)) void erasure_code_process(
    uint8_t *table,
    int input_count,
    int output_count,
    const uint8_t *restrict*restrict input_shard_data,
    uint8_t *restrict*restrict output_shard_data,
    size_t shard_length
) {
    assert(shard_length % 32 == 0);
    const uint64_t *restrict*restrict input64 = (const uint64_t *restrict*restrict)input_shard_data;
    uint64_t *restrict*restrict output64 = (uint64_t *restrict*restrict)output_shard_data;
    size_t segments = shard_length / 32;
    //const uint64_t *input64_i = input64[i];
    for (int j = 0; j < output_count; j++) {
    //uint64_t *output64_j = output64[j];
        for (int segment = 0; segment < segments; segment++) {
            uint64_t *restrict output_ptr = output64[j] + segment * 4;
            uint64_t a = 0, b = 0, c = 0, d = 0;
            for (int i = 0; i < input_count; i++) {
                const uint64_t *restrict input_ptr = input64[i] + segment * 4;
                for (int bit_in = 0; bit_in < 4; bit_in++) {
                    uint8_t mask = table[i * output_count * 4 + j * 4 + bit_in];
                    uint64_t bit = input_ptr[bit_in];
                    //a ^= bit * (mask & 1);
                    //b ^= bit * ((mask >> 1) & 1);
                    //c ^= bit * ((mask >> 2) & 1);
                    //d ^= bit * ((mask >> 3) & 1);
                    a ^= (mask & 1) ? bit : 0;
                    b ^= (mask & 2) ? bit : 0;
                    c ^= (mask & 4) ? bit : 0;
                    d ^= (mask & 8) ? bit : 0;
                    /*
                    if (mask & 1) output_ptr[0] ^= input_ptr[bit_in];
                    if (mask & 2) output_ptr[1] ^= input_ptr[bit_in];
                    if (mask & 4) output_ptr[2] ^= input_ptr[bit_in];
                    if (mask & 8) output_ptr[3] ^= input_ptr[bit_in];
                    */
                }
            }
            output_ptr[0] = a;
            output_ptr[1] = b;
            output_ptr[2] = c;
            output_ptr[3] = d;
        }
    }
}

void erasure_code_process_bad(
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


#define BYTES (2 * 1024 * 1024)

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
