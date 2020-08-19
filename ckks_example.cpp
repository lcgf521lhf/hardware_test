#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#define Enc_ntt_size 3
#define Coeff_modulus_size 4
#define Coeff_count 8192

using namespace std;

typedef unsigned long long u64;

inline unsigned char add_uint64(u64 operand1, u64 operand2, u64 *result)
{
    *result = operand1 + operand2;
    return static_cast<unsigned char>(*result < operand1);
}

inline void multiply_uint64(u64 operand1, u64 operand2, u64 *result128)
{
    auto operand1_coeff_right = operand1 & 0x00000000FFFFFFFFULL;
    auto operand2_coeff_right = operand2 & 0x00000000FFFFFFFFULL;
    operand1 >>= 32;
    operand2 >>= 32;

    auto middle1 = operand1 * operand2_coeff_right;
    u64 middle;
    auto left =
        operand1 * operand2 + (static_cast<u64>(add_uint64(middle1, operand2 * operand1_coeff_right, &middle)) << 32);
    auto right = operand1_coeff_right * operand2_coeff_right;
    auto temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFFULL);

    result128[1] = static_cast<u64>(left + (middle >> 32) + (temp_sum >> 32));
    result128[0] = static_cast<u64>((temp_sum << 32) | (right & 0x00000000FFFFFFFFULL));
}

inline void multiply_uint64_hw64(u64 operand1, u64 operand2, u64 *hw64)
{
    auto operand1_coeff_right = operand1 & 0x00000000FFFFFFFFULL;
    auto operand2_coeff_right = operand2 & 0x00000000FFFFFFFFULL;
    operand1 >>= 32;
    operand2 >>= 32;

    auto middle1 = operand1 * operand2_coeff_right;
    u64 middle;
    auto left =
        operand1 * operand2 + (static_cast<u64>(add_uint64(middle1, operand2 * operand1_coeff_right, &middle)) << 32);
    auto right = operand1_coeff_right * operand2_coeff_right;
    auto temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFFULL);

    *hw64 = static_cast<u64>(left + (middle >> 32) + (temp_sum >> 32));
}

int main()
{
    u64 *plain_ntt = (u64 *)malloc(Coeff_modulus_size * Coeff_count * sizeof(u64));
    if (plain_ntt == nullptr)
    {
        return -1;
    }

    ifstream in1;
    in1.open("plain.txt", ios::in);
    for (size_t i = 0; i < 3 * Coeff_count; i++)
    {
        in1 >> plain_ntt[i];
    }

    u64 *encrypted_ntt = (u64 *)malloc(Enc_ntt_size * Coeff_modulus_size * Coeff_count * sizeof(u64));
    if (encrypted_ntt == nullptr)
    {
        return -1;
    }

    ifstream in2;
    in2.open("enc.txt", ios::in);
    for (size_t i = 0; i < 2 * 3 * Coeff_count; i++)
    {
        in2 >> encrypted_ntt[i];
    }

    u64 *mod = (u64 *)malloc(9 * sizeof(u64));
    ifstream in3;
    in3.open("modulus.txt", ios::in);
    for (size_t i = 0; i < 9; i++)
    {
        in3 >> mod[i];
    }

    auto ts = std::chrono::high_resolution_clock::now();
    int coeff_modulus_size = 3;
    int coeff_count = 8192;
    for (size_t i = 0; i < Enc_ntt_size; i++)
    {
        if (i >= 2)
        {
            break;
        }
        for (size_t j = 0; j < Coeff_modulus_size; j++)
        {
            if (j >= coeff_modulus_size)
            {
                break;
            }
            const u64 modulus_value = mod[3 * j];
            const u64 const_ratio_0 = mod[3 * j + 1];
            const u64 const_ratio_1 = mod[3 * j + 2];
            for (size_t k = 0; k < Coeff_count; k++)
            {
                if (k >= coeff_count)
                {
                    break;
                }
                // Reduces z using base 2^64 Barrett reduction
                u64 z[2], tmp1, tmp2[2], tmp3, carry;
                multiply_uint64(
                    encrypted_ntt[i * coeff_modulus_size * coeff_count + j * coeff_count + k],
                    plain_ntt[j * coeff_count + k], z);

                // Multiply input and const_ratio
                // Round 1
                multiply_uint64_hw64(z[0], const_ratio_0, &carry);
                multiply_uint64(z[0], const_ratio_1, tmp2);
                tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, &tmp1);

                // Round 2
                multiply_uint64(z[1], const_ratio_0, tmp2);
                carry = tmp2[1] + add_uint64(tmp1, tmp2[0], &tmp1);

                // This is all we care about
                tmp1 = z[1] * const_ratio_1 + tmp3 + carry;

                // Barrett subtraction
                tmp3 = z[0] - tmp1 * modulus_value;

                // Claim: One more subtraction is enough
                encrypted_ntt[i * coeff_modulus_size * coeff_count + j * coeff_count + k] =
                    tmp3 - (modulus_value & static_cast<u64>(-static_cast<signed long long>(tmp3 >= modulus_value)));
            }
        }
    }

    auto te = std::chrono::high_resolution_clock::now();
    cout << "[Serial] " << std::chrono::duration_cast<std::chrono::nanoseconds>(te - ts).count() << " ns - \n";

    ofstream out3;
    out3.open("res2.txt", ios::in | ios::out);
    for (size_t i = 0; i < 2 * coeff_modulus_size * coeff_count; i++)
    {
        out3 << encrypted_ntt[i] << endl;
    }
}