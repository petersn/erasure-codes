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
__attribute__((alwaysinline)) inline __m256i cond_xor(
    __m256i accum,
    uint32_t masks,
    int bit_pos,
    __m256i bit
) {
    if (bit_pos > 7) {
        masks >>= bit_pos - 7;
    } else {
        masks <<= 7 - bit_pos;
    }
    __m256i full_mask = _mm256_set1_epi8(masks);
    __m256i xor = _mm256_xor_si256(accum, bit);
    return _mm256_blendv_epi8(accum, xor, full_mask);
    //return _mm256_xor_si256(accum, (masks & (1u << bit_pos)) ? bit : _mm256_setzero_si256());
}

__attribute__((noinline)) void erasure_code_process(
    uint8_t *restrict table,
    int input_count,
    int output_count,
    const uint8_t *restrict*restrict input_shard_data,
    uint8_t *restrict*restrict output_shard_data,
    size_t shard_length
) {
    // Check the alignment on the input and output shards.
    for (int i = 0; i < input_count; i++)
        assert(((uintptr_t)input_shard_data[i] & 0x1f) == 0);
    for (int i = 0; i < output_count; i++)
        assert(((uintptr_t)output_shard_data[i] & 0x1f) == 0);
    assert(shard_length % 32 == 0);
    const __m256i *restrict*restrict input64 = (const __m256i *restrict*restrict)input_shard_data;
    __m256i *restrict*restrict output64 = (__m256i *restrict*restrict)output_shard_data;
    uint32_t *restrict table32 = (uint32_t *restrict)table;
    size_t segments = shard_length / 128;
    //const __m256i *input64_i = input64[i];
    for (int j = 0; j < output_count; j++) {
    //__m256i *output64_j = output64[j];
        for (int segment = 0; segment < segments; segment++) {
            __m256i *restrict output_ptr = output64[j] + segment * 4;
            __m256i a = _mm256_setzero_si256(), b = _mm256_setzero_si256(), c = _mm256_setzero_si256(), d = _mm256_setzero_si256();
            for (int i = 0; i < input_count; i++) {
                const __m256i *restrict input_ptr = input64[i] + segment * 4;
                uint32_t masks = table32[i * output_count * 4 + j * 4];
                __m256i bit0 = input_ptr[0];
                __m256i bit1 = input_ptr[1];
                __m256i bit2 = input_ptr[2];
                __m256i bit3 = input_ptr[3];
                //a = _mm256_xor_si256(a, (masks & (1u <<  0)) ? bit0 : _mm256_setzero_si256());
                //b = _mm256_xor_si256(b, (masks & (1u <<  1)) ? bit0 : _mm256_setzero_si256());
                //c = _mm256_xor_si256(c, (masks & (1u <<  2)) ? bit0 : _mm256_setzero_si256());
                //d = _mm256_xor_si256(d, (masks & (1u <<  3)) ? bit0 : _mm256_setzero_si256());
                //a = _mm256_xor_si256(a, (masks & (1u <<  4)) ? bit1 : _mm256_setzero_si256());
                //b = _mm256_xor_si256(b, (masks & (1u <<  5)) ? bit1 : _mm256_setzero_si256());
                //c = _mm256_xor_si256(c, (masks & (1u <<  6)) ? bit1 : _mm256_setzero_si256());
                //d = _mm256_xor_si256(d, (masks & (1u <<  7)) ? bit1 : _mm256_setzero_si256());
                //a = _mm256_xor_si256(a, (masks & (1u <<  8)) ? bit2 : _mm256_setzero_si256());
                //b = _mm256_xor_si256(b, (masks & (1u <<  9)) ? bit2 : _mm256_setzero_si256());
                //c = _mm256_xor_si256(c, (masks & (1u << 10)) ? bit2 : _mm256_setzero_si256());
                //d = _mm256_xor_si256(d, (masks & (1u << 11)) ? bit2 : _mm256_setzero_si256());
                //a = _mm256_xor_si256(a, (masks & (1u << 12)) ? bit3 : _mm256_setzero_si256());
                //b = _mm256_xor_si256(b, (masks & (1u << 13)) ? bit3 : _mm256_setzero_si256());
                //c = _mm256_xor_si256(c, (masks & (1u << 14)) ? bit3 : _mm256_setzero_si256());
                //d = _mm256_xor_si256(d, (masks & (1u << 15)) ? bit3 : _mm256_setzero_si256());
                a = cond_xor(a, masks,  0, bit0);
                b = cond_xor(b, masks,  1, bit0);
                c = cond_xor(c, masks,  2, bit0);
                d = cond_xor(d, masks,  3, bit0);
                a = cond_xor(a, masks,  4, bit1);
                b = cond_xor(b, masks,  5, bit1);
                c = cond_xor(c, masks,  6, bit1);
                d = cond_xor(d, masks,  7, bit1);
                a = cond_xor(a, masks,  8, bit2);
                b = cond_xor(b, masks,  9, bit2);
                c = cond_xor(c, masks, 10, bit2);
                d = cond_xor(d, masks, 11, bit2);
                a = cond_xor(a, masks, 12, bit3);
                b = cond_xor(b, masks, 13, bit3);
                c = cond_xor(c, masks, 14, bit3);
                d = cond_xor(d, masks, 15, bit3);
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
        out[0][3*i + 0] = data[0][3*i + 0] ^ data[0][3*i + 1] ^ data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[1][3*i + 2] ^ data[2][3*i + 1] ^ data[2][3*i + 2] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 0];
        out[0][3*i + 1] = data[0][3*i + 0] ^ data[1][3*i + 0] ^ data[2][3*i + 0] ^ data[2][3*i + 1] ^ data[3][3*i + 0] ^ data[3][3*i + 1] ^ data[4][3*i + 1];
        out[0][3*i + 2] = data[0][3*i + 0] ^ data[0][3*i + 1] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[2][3*i + 0] ^ data[2][3*i + 1] ^ data[2][3*i + 2] ^ data[3][3*i + 0] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 2];

        out[1][3*i + 0] = data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[1][3*i + 2] ^ data[2][3*i + 2] ^ data[3][3*i + 0] ^ data[3][3*i + 2] ^ data[4][3*i + 0];
        out[1][3*i + 1] = data[0][3*i + 0] ^ data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[1][3*i + 2] ^ data[2][3*i + 0] ^ data[2][3*i + 2] ^ data[3][3*i + 0] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 1];
        out[1][3*i + 2] = data[0][3*i + 1] ^ data[1][3*i + 1] ^ data[1][3*i + 2] ^ data[2][3*i + 1] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 2];

        out[2][3*i + 0] = data[0][3*i + 1] ^ data[1][3*i + 0] ^ data[1][3*i + 1] ^ data[2][3*i + 0] ^ data[2][3*i + 1] ^ data[3][3*i + 1] ^ data[4][3*i + 0];
        out[2][3*i + 1] = data[0][3*i + 1] ^ data[0][3*i + 2] ^ data[1][3*i + 2] ^ data[2][3*i + 2] ^ data[3][3*i + 1] ^ data[3][3*i + 2] ^ data[4][3*i + 1];
        out[2][3*i + 2] = data[0][3*i + 0] ^ data[0][3*i + 2] ^ data[1][3*i + 0] ^ data[2][3*i + 0] ^ data[3][3*i + 0] ^ data[3][3*i + 2] ^ data[4][3*i + 2];
    }
}

void* aligned_malloc(size_t size, size_t alignment) {
    void* ptr = NULL;
    if (posix_memalign(&ptr, alignment, size) != 0) {
        return NULL;
    }
    return ptr;
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
        shards[i] = aligned_malloc(BYTES, 32);
        for (int j = 0; j < BYTES; j++) {
            shards[i][j] = i * BYTES + j;
        }
    }

    for (int j = 0; j < 1000; j++) {
        erasure_code_process(table, 5, 3, (const uint8_t**)shards, shards + 5, BYTES);
    }
}
