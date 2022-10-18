#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

const uint8_t exp_table[255 * 2] = {
    1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19, 53, 95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34, 102, 170, 229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112, 144, 171, 230, 49, 83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104, 184, 211, 110, 178, 205, 76, 212, 103, 169, 224, 59, 77, 215, 98, 166, 241, 8, 24, 40, 120, 136, 131, 158, 185, 208, 107, 189, 220, 127, 129, 152, 179, 206, 73, 219, 118, 154, 181, 196, 87, 249, 16, 48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163, 254, 25, 43, 125, 135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160, 251, 22, 58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65, 195, 94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218, 117, 159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223, 122, 142, 137, 128, 155, 182, 193, 88, 232, 35, 101, 175, 234, 37, 111, 177, 200, 67, 197, 84, 252, 31, 33, 99, 165, 244, 7, 9, 27, 45, 119, 153, 176, 203, 70, 202, 69, 207, 74, 222, 121, 139, 134, 145, 168, 227, 62, 66, 198, 81, 243, 14, 18, 54, 90, 238, 41, 123, 141, 140, 143, 138, 133, 148, 167, 242, 13, 23, 57, 75, 221, 124, 132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246,
    1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19, 53, 95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34, 102, 170, 229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112, 144, 171, 230, 49, 83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104, 184, 211, 110, 178, 205, 76, 212, 103, 169, 224, 59, 77, 215, 98, 166, 241, 8, 24, 40, 120, 136, 131, 158, 185, 208, 107, 189, 220, 127, 129, 152, 179, 206, 73, 219, 118, 154, 181, 196, 87, 249, 16, 48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163, 254, 25, 43, 125, 135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160, 251, 22, 58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65, 195, 94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218, 117, 159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223, 122, 142, 137, 128, 155, 182, 193, 88, 232, 35, 101, 175, 234, 37, 111, 177, 200, 67, 197, 84, 252, 31, 33, 99, 165, 244, 7, 9, 27, 45, 119, 153, 176, 203, 70, 202, 69, 207, 74, 222, 121, 139, 134, 145, 168, 227, 62, 66, 198, 81, 243, 14, 18, 54, 90, 238, 41, 123, 141, 140, 143, 138, 133, 148, 167, 242, 13, 23, 57, 75, 221, 124, 132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246
};

const uint8_t log_table[256] = {
    0, 255, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223, 3, 100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105, 28, 193, 125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114, 154, 201, 9, 120, 101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130, 69, 53, 147, 218, 142, 150, 143, 219, 189, 54, 208, 206, 148, 19, 92, 210, 241, 64, 70, 131, 56, 102, 221, 253, 48, 191, 6, 139, 98, 179, 37, 226, 152, 34, 136, 145, 16, 126, 110, 72, 195, 163, 182, 30, 66, 58, 107, 40, 84, 250, 133, 61, 186, 43, 121, 10, 21, 155, 159, 94, 202, 78, 212, 172, 229, 243, 115, 167, 87, 175, 88, 168, 80, 244, 234, 214, 116, 79, 174, 233, 213, 231, 230, 173, 232, 44, 215, 117, 122, 235, 22, 11, 245, 89, 203, 95, 176, 156, 169, 81, 160, 127, 12, 246, 111, 23, 196, 73, 236, 216, 67, 31, 45, 164, 118, 123, 183, 204, 187, 62, 90, 251, 96, 177, 134, 59, 82, 161, 108, 170, 85, 41, 157, 151, 178, 135, 144, 97, 190, 220, 252, 188, 149, 207, 205, 55, 63, 91, 209, 83, 57, 132, 60, 65, 162, 109, 71, 20, 42, 158, 93, 86, 242, 211, 171, 68, 17, 146, 217, 35, 32, 46, 137, 180, 124, 184, 38, 119, 153, 227, 165, 103, 74, 237, 222, 197, 49, 254, 24, 13, 99, 140, 128, 192, 247, 112, 7
};

//#define LOOKUP_TABLES
//#define SHORT_EXP_TABLE

#ifdef LOOKUP_TABLES
uint8_t field_mul(uint8_t a, uint8_t b) {
    if (a == 0 || b == 0)
        return 0;
#ifdef SHORT_EXP_TABLE
    return exp_table[(((int)log_table[a]) + ((int)log_table[b])) % 255];
#else
    return exp_table[((int)log_table[a]) + ((int)log_table[b])];
#endif
}
#else
uint8_t field_mul(uint8_t a, uint8_t b) {
    uint16_t x = a;
    uint8_t result = 0;
    for (int i = 0; i < 8; i++) {
        if (b & 1)
            result ^= x;
        x <<= 1;
        if (x & 0x100)
            x ^= 0x1b;
        b >>= 1;
    }
    return result;
}
#endif

uint8_t field_inv(uint8_t a) {
    // Raise to the 254th power
    uint16_t accum = 1;
    for (int i = 0; i < 7; i++) {
        a = field_mul(a, a);
        accum = field_mul(accum, a);
    }
    return accum;
}

void erasure_code_build_table(
    uint8_t *table,
    int input_count,
    int output_count,
    const int *input_x,
    const int *output_x
) {
    assert(input_count + output_count < 256);
    for (int i = 0; i < input_count; i++) {
        for (int j = 0; j < output_count; j++) {
            uint8_t accum = 1;
            for (int p = 0; p < input_count; p++) {
                if (p == i)
                    continue;
                accum = field_mul(accum, field_mul(
                    output_x[j] ^ input_x[p],
                    field_inv(input_x[i] ^ input_x[p])
                ));
            }
            table[i * output_count + j] = accum;
        }
    }
}

__attribute__((noinline)) void erasure_code_process(
    uint8_t *table,
    int input_count,
    int output_count,
    const uint8_t **input_shard_data,
    uint8_t **output_shard_data,
    size_t shard_length
) {
    for (int i = 0; i < output_count; i++) {
        for (size_t j = 0; j < shard_length; j++) {
            output_shard_data[i][j] = 0;
        }
    }
    for (int i = 0; i < input_count; i++) {
        for (int j = 0; j < output_count; j++) {
            uint8_t coeff = table[j + i * output_count];
            const uint8_t *input_ptr = input_shard_data[i];
            uint8_t *output_ptr = output_shard_data[j];
            for (size_t k = 0; k < shard_length; k++)
                output_ptr[k] ^= field_mul(coeff, input_ptr[k]);
        }
    }
}

#define BYTES (1024 * 1024)

int main() {
    // Configure to encode 5 data shards into 3 parity shards.
    uint8_t table[5 * 3];
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

/*
int main() {
    erasure_code_context_t ctx;
    erasure_code_init(&ctx, 5, 3);

    // Configure to encode 5 data shards into 3 parity shards.
    int input_x[5] = {0, 1, 2, 3, 4};
    int output_x[3] = {5, 6, 7};
    erasure_code_configure(&ctx, input_x, output_x);

    uint8_t *shards[8 + 3];
    for (int i = 0; i < 8 + 3; i++) {
        shards[i] = malloc(16);
        for (int j = 0; j < 16; j++)
            shards[i][j] = i + j;
    }
    erasure_code_process(&ctx, shards, shards + 5, 16);

    // Print out the results.
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 16; j++) {
            printf("%02x", shards[i][j]);
        }
        printf("\n");
    }

    // Configure to recover from the last five shards.
    int input_x2[5] = {3, 4, 5, 6, 7};
    int output_x2[3] = {0, 1, 2};
    erasure_code_configure(&ctx, input_x2, output_x2);
    erasure_code_process(&ctx, shards + 3, shards + 8, 16);

    printf("Reconstructed first three shards:\n");
    for (int i = 8; i < 8 + 3; i++) {
        for (int j = 0; j < 16; j++) {
            printf("%02x", shards[i][j]);
        }
        printf("\n");
    }

    erasure_code_free(&ctx);
    return 0;
}
*/
