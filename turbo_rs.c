#include <stdint.h>
#include <unistd.h>
#include <assert.h>

#if defined(_M_X64) && !defined(__x86_64__)
#define __x86_64__
#endif

#if defined(_M_ARM64) && !defined(__aarch64__)
#define __aarch64__
#endif

#ifndef NO_JIT
#include <sys/mman.h>
#endif

#if defined(__x86_64__) && !defined(NO_AVX2)
#define REED_SOLOMON_AVX2
#include <x86intrin.h>
#endif

#if defined(__aarch64__) && !defined(NO_NEON)
#define REED_SOLOMON_NEON
// TODO
#endif

static inline int get_degree(uint64_t poly) {
    return 63 - __builtin_clzll(poly);
}

static inline uint32_t field_mul(uint64_t poly, uint32_t a, uint32_t b) {
    uint64_t x = a;
    uint32_t result = 0;
    for (int i = 0; i < get_degree(poly); i++) {
        if (b & 1)
            result ^= x;
        x <<= 1;
        if (x & (1ull << get_degree(poly)))
            x ^= poly;
        b >>= 1;
    }
    return result;
}

static inline uint32_t field_inv(uint64_t poly, uint32_t a) {
    // Raise to the (2**d - 1)th power
    uint16_t accum = 1;
    for (int i = 0; i < get_degree(poly) - 1; i++) {
        a = field_mul(poly, a, a);
        accum = field_mul(poly, accum, a);
    }
    return accum;
}

// Takes two polynomials over GF(2) stored little-endian with one bit per byte.
// Each len must be the index of the highest bit set + 1.
// Returns a pointer to the gcd polynomial, and sets *len_out to its length.
static uint8_t poly_gcd(
    int len1,
    uint8_t *poly1,
    int len2,
    uint8_t *poly2,
    int *len_out
) {
    while (len1 && len2) {
        if (len1 > len2) {
            // Swap the two polynomials.
            int tmp_len = len1;
            len1 = len2;
            len2 = tmp_len;
            uint8_t *tmp_poly = poly1;
            poly1 = poly2;
            poly2 = tmp_poly;
        }
        int shift_amount = len2 - len1;
        len2 = 0;
        for (int i = 0; i < len1; i++) {
            poly2[i + shift_amount] ^= poly1[i];
            if (poly2[i + shift_amount])
                len2 = i + shift_amount + 1;
        }
    }
    *len_out = len1 ? len1 : len2;
    return len1 ? poly1 : poly2;
}

/*
bool check_is_irreducible(uint64_t poly) {
    int degree = get_degree(poly);
    // First check if poly divides x**(2**deg) + x.
    uint8_t *test_poly = (uint8_t *)malloc(degree + 1);
    for (int i = 0; i <= degree; i++)
        test_poly[i] = (poly >> i) & 1;
    uint8_t *big_poly = (uint8_t *)malloc((1 << degree) + 1);
    big_poly[1] = 1;
    int big_poly_len;
    poly_gcd(degree + 1, test_poly, 1ull << get_degree(poly) + 1, big_poly, &big_poly_len);
}
*/

int check_is_irreducible(uint64_t poly) {
    // TODO: Implement this.
    return 1;
}

typedef struct reed_solomon_t {
    uint64_t polynomial;
    int input_count;
    int output_count;
    uint32_t *table;
    uint32_t *table2;
    void *jit_code;
    size_t jit_code_length;
} reed_solomon_t;

typedef enum rs_processing_mode_t {
    RS_PROCESSING_MODE_FASTEST_AVAILABLE,
    RS_PROCESSING_MODE_BASIC,
    RS_PROCESSING_MODE_VECTORIZED,
    RS_PROCESSING_MODE_JIT,
    RS_PROCESSING_MODE_JIT_VECTORIZED,
} rs_processing_mode_t;

int alignment_requirement(
    uint64_t polynomial,
    rs_processing_mode_t mode
) {
    switch (mode) {
        default:
    }
    assert(0);
}

int byte_multiple_requirement(
    uint64_t polynomial,
    rs_processing_mode_t mode
) {
    switch (mode) {
        default:
    }
    assert(0);
}

void reed_solomon_init(
    reed_solomon_t *ctx,
    uint64_t polynomial,
    int input_count,
    int output_count
) {
    int degree = get_degree(polynomial);
    assert(degree <= 32);
    assert(input_count + output_count <= (1 << degree));
    assert(check_is_irreducible(polynomial));

    ctx->polynomial = polynomial;
    ctx->input_count = input_count;
    ctx->output_count = output_count;
    ctx->table = (uint32_t*) malloc(
        sizeof(uint32_t) * input_count * output_count * degree
    );
    ctx->table2 = (uint32_t*) malloc(
        sizeof(uint32_t) * input_count * output_count * degree
    );
    ctx->jit_code = NULL;
    ctx->jit_code_length = 0;
}

void reed_solomon_free(reed_solomon_t *ctx) {
    free(ctx->table);
    ctx->table = NULL;
    if (ctx->jit_code != NULL)
        munmap(ctx->jit_code, ctx->jit_code_length);
    ctx->jit_code = NULL;
}

void reed_solomon_configure(
    reed_solomon_t *ctx,
    const int *input_x,
    const int *output_x
) {
    int degree = get_degree(ctx->polynomial);
    for (int i = 0; i < ctx->input_count; i++) {
        assert(0 <= input_x[i] && input_x[i] < (1 << degree));
        for (int j = 0; j < ctx->output_count; j++) {
            assert(0 <= output_x[j] && output_x[j] < (1 << degree));
            uint32_t accum = 1;
            for (int p = 0; p < ctx->input_count; p++) {
                if (p == i)
                    continue;
                accum = field_mul(ctx->polynomial, accum, field_mul(
                    ctx->polynomial,
                    output_x[j] ^ input_x[p],
                    field_inv(ctx->polynomial, input_x[i] ^ input_x[p])
                ));
            }
            ctx->table2[i * ctx->output_count + j] = accum;
            for (int input_bit = 0; input_bit < degree; input_bit++) {
                ctx->table[i * ctx->output_count * degree + j * degree + input_bit] =
                    field_mul(ctx->polynomial, accum, 1ull << input_bit);
            }
        }
    }
    // JIT compile the code.
    if (ctx->jit_code != NULL)
        munmap(ctx->jit_code, ctx->jit_code_length);
    ctx->jit_code = mmap(
        NULL,
        4096,
        PROT_READ | PROT_WRITE | PROT_EXEC, MAP_PRIVATE | MAP_ANONYMOUS,
        -1,
        0
    );
}

void reed_solomon_process_jit(
) {
    
}

void reed_solomon_process_no_jit(
    reed_solomon_t *ctx,
    const char * restrict * input_shard_data,
    char * restrict * output_shard_data,
    size_t shard_length
) {
}

void reed_solomon_process_no_jit_bytewise(
    reed_solomon_t *ctx,
    const char * restrict * input_shard_data,
    char * restrict * output_shard_data,
    size_t shard_length
) {
    for (int i = 0; i < ctx->output_count; i++) {
        for (size_t j = 0; j < shard_length; j++) {
            output_shard_data[i][j] = 0;
        }
    }
    for (int i = 0; i < ctx->input_count; i++) {
        for (int j = 0; j < ctx->output_count; j++) {
            uint8_t coeff = ctx->table2[i * ctx->output_count + j];
            const uint8_t *input_ptr = input_shard_data[i];
            uint8_t *output_ptr = output_shard_data[j];
            for (size_t k = 0; k < shard_length; k++)
                output_ptr[k] ^= field_mul(ctx->polynomial, coeff, input_ptr[k]);
        }
    }
}

#ifdef REED_SOLOMON_AVX2
void reed_solomon_process_no_jit_avx2(
    reed_solomon_t *ctx,
    const char * restrict * input_shard_data,
    char * restrict * output_shard_data,
    size_t shard_length
) {
    int degree = get_degree(ctx->polynomial);
    // Check the alignment on the input and output shards.
    for (int i = 0; i < ctx->input_count; i++)
        assert(((uintptr_t) input_shard_data[i] & 0x1f) == 0);
    for (int i = 0; i < ctx->output_count; i++)
        assert(((uintptr_t) output_shard_data[i] & 0x1f) == 0);
    assert(shard_length % (32 * degree) == 0);
    size_t segments = shard_length / (32 * degree);

    int input_count = ctx->input_count, output_count = ctx->output_count;
    const __m256i **input_vec = (const __m256i **) input_shard_data;
    __m256i ** output_vec = (__m256i **) output_shard_data;
    uint32_t *table = ctx->table;
    __m256i scratch[degree];
    for (int j = 0; j < output_count; j++) {
        for (size_t segment = 0; segment < segments; segment++) {
            __m256i * restrict output_ptr = output_vec[j] + segment * degree;
            for (int b = 0; b < degree; b++)
                scratch[b] = _mm256_setzero_si256();
            for (int i = 0; i < input_count; i++) {
                const __m256i *restrict input_ptr = input_vec[i] + segment * degree;
                for (int input_bit = 0; input_bit < degree; input_bit++) {
                    uint32_t mask = table[i * output_count * degree + j * degree + input_bit];
                    for (int k = 0; k < degree; k++) {
                        uint32_t mask_bit = ((mask >> k) & 1) << 7;
                        __m256i full_mask = _mm256_set1_epi8(mask_bit);
                        __m256i xor = _mm256_xor_si256(scratch[k], input_ptr[input_bit]);
                        scratch[k] = _mm256_blendv_epi8(scratch[k], xor, full_mask);
                    }
                }
            }
            for (int b = 0; b < degree; b++)
                output_ptr[b] = scratch[b];
        }
    }
}
#endif

void* aligned_malloc(size_t size, size_t alignment) {
    void* ptr = NULL;
    if (posix_memalign(&ptr, alignment, size) != 0)
        return NULL;
    return ptr;
}

#include <stdio.h>

//#define BYTES (1024 * 1024)
#define BYTES 96

int main() {
    // Example usage.
    reed_solomon_t ctx;
    reed_solomon_init(&ctx, 0b1011/*0x11d*/, 5, 3);
    //reed_solomon_init(&ctx, 0x11b, 5, 3);
    int input_x[5] = {0, 1, 2, 3, 4};
    int output_x[3] = {5, 6, 7};
    reed_solomon_configure(&ctx, input_x, output_x);

    // Print table2.
    for (int i = 0; i < ctx.input_count; i++) {
        for (int j = 0; j < ctx.output_count; j++) {
            printf("%02x ", ctx.table2[i * ctx.output_count + j]);
        }
        printf("\n");
    }
    printf("DONE\n");

    char *shards[8];
    for (int i = 0; i < 8; i++) {
        shards[i] = aligned_malloc(BYTES, 32);
        for (int j = 0; j < BYTES; j++) {
            shards[i][j] = i < 5 ? i + j : 0;
        }
    }
    reed_solomon_process_no_jit_avx2(
        &ctx, (const char**)shards, shards + 5, BYTES
    );
    // Print out the first few bytes of the shards.
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 16; j++) {
            printf("%02x ", (unsigned char) shards[i][j]);
        }
        printf("\n");
    }
    //return 0;

    // Reconfigure to recover.
    int input_x2[5] = {0, 2, 4, 5, 7};
    int output_x2[3] = {1, 3, 6};
    reed_solomon_configure(&ctx, input_x2, output_x2);
    // Print table2.
    for (int i = 0; i < ctx.input_count; i++) {
        for (int j = 0; j < ctx.output_count; j++) {
            printf("%02x ", ctx.table2[i * ctx.output_count + j]);
        }
        printf("\n");
    }
    printf("DONE\n");
    char* shards2[8];
    for (int i = 0; i < 5; i++) {
        shards2[i] = shards[input_x2[i]];
    }
    for (int i = 5; i < 8; i++) {
        shards2[i] = aligned_malloc(BYTES, 32);
        memset(shards2[i], 0, BYTES);
    }

    reed_solomon_process_no_jit_avx2(
        &ctx, (const char**)shards2, shards2 + 5, BYTES
    );
    // Compare the recovered shards.
    int bad = 0;
    for (int i = 0; i < 3; i++) {
        printf("======\nOriginal:  ");
        // Compare shards[output_x2[i]] and shards2[5 + i]
        for (int j = 0; j < 16; j++) {
            printf("%02x ", (unsigned char) shards[output_x2[i]][j]);
            bad |= shards[output_x2[i]][j] != shards2[5 + i][j];
        }
        printf("\nRecovered: ");
        for (int j = 0; j < 16; j++) {
            printf("%02x ", (unsigned char) shards2[5 + i][j]);
        }
        printf("\n======\n");
    }
    printf("BAD: %d\n", bad);

    //reed_solomon_process_jit(&ctx, RS_PROCESSING_MODE_AVX2);

}
