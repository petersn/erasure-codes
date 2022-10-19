#include <stdint.h>
#include <unistd.h>
#include <assert.h>

#if defined(__x86_64__) || defined(_M_X64)
#define __x86_64__
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
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

typedef struct reed_solomon_t {
    uint64_t polynomial;
    int input_count;
    int output_count;
    uint32_t *table;
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
        case RS_PROCESSING_MODE_FASTEST_AVAILABLE:
        case RS_PROCESSING_MODE_BASIC:
            return 1;
        case RS_PROCESSING_MODE_VECTORIZED:
        case RS_PROCESSING_MODE_JIT:
        case RS_PROCESSING_MODE_JIT_VECTORIZED:
            return 4;
    }
    assert(0);
}

int 

void reed_solomon_init(
    reed_solomon_t *ctx,
    uint64_t polynomial,
    int input_count,
    int output_count
) {
    int degree = get_degree(polynomial);
    assert(degree <= 32);
    assert(input_count + output_count <= (1 << degree));

    ctx->polynomial = polynomial;
    ctx->input_count = input_count;
    ctx->output_count = output_count;
    ctx->table = malloc(sizeof(uint32_t) * input_count * output_count * degree);
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
        assert(0 <= input_x[i] && input_x[i] < degree);
        for (int j = 0; j < ctx->output_count; j++) {
            assert(0 <= output_x[j] && output_x[j] < degree);
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
            for (int bit = 0; bit < degree; bit++) {
                ctx->table[i * ctx->output_count * degree + j * degree + bit] =
                    field_mul(ctx->polynomial, accum, 1ull << bit);
            }
        }
    }
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
        for (int segment = 0; segment < segments; segment++) {
            __m256i * restrict output_ptr = output_vec[j] + segment * degree;
            for (int i = 0; i < degree; i++)
                scratch[i] = _mm256_setzero_si256();
            for (int i = 0; i < input_count; i++) {
                const __m256i *restrict input_ptr = input_vec[i] + segment * degree;
                uint32_t mask = table[i * output_count * 4 + j * 4];
                for (int k = 0; k < degree; k++) {
                    uint32_t bit = ((mask >> k) & 1) << 7;
                    __m256i full_mask = _mm256_set1_epi8(bit);
                    __m256i xor = _mm256_xor_si256(scratch[k], input_ptr[k]);
                    scratch[k] = _mm256_blendv_epi8(scratch[k], xor, full_mask);
                }
            }
            for (int i = 0; i < degree; i++)
                output_ptr[i] = scratch[i];
        }
    }
}
#endif

void reed_solomon_process(
    reed_solomon_t *ctx,
    const char * restrict * input_shard_data,
    char * restrict * output_shard_data,
    size_t shard_length

//mmap(NULL, 4096, PROT_READ | PROT_WRITE | PROT_EXEC, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

int main() {
    
}
