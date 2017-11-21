#include <iostream>

//using namespace std;
#include <cryptoTools/Common/Defines.h>
using namespace osuCrypto;

#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Endpoint.h>
#include <numeric>
#include <cryptoTools/Common/Log.h>
#include <cryptoTools/Common/Timer.h>

#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtReceiver.h"
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtSender.h"
#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"

#include <cryptoTools/gsl/span>
#include <cryptoTools/Common/ByteStream.h>
#include <cryptoTools/Common/Matrix.h>

using namespace std;
static void sse_trans(uint8_t *out, uint8_t const *inp, int nrows, int ncols) {
#   define INP(x,y) inp[(x)*ncols/8 + (y)/8]
#   define OUT(x,y) out[(y)*nrows/8 + (x)/8]
	int rr, cc, i, h;
	union { __m128i x; uint8_t b[16]; } tmp;
	__m128i vec;
	assert(nrows % 8 == 0 && ncols % 8 == 0);

	// Do the main body in 16x8 blocks:
	for (rr = 0; rr <= nrows - 16; rr += 16) {
		for (cc = 0; cc < ncols; cc += 8) {
			vec = _mm_set_epi8(
				INP(rr + 15, cc), INP(rr + 14, cc), INP(rr + 13, cc), INP(rr + 12, cc), INP(rr + 11, cc), INP(rr + 10, cc), INP(rr + 9, cc),
				INP(rr + 8, cc), INP(rr + 7, cc), INP(rr + 6, cc), INP(rr + 5, cc), INP(rr + 4, cc), INP(rr + 3, cc), INP(rr + 2, cc), INP(rr + 1, cc),
				INP(rr + 0, cc));
			for (i = 8; --i >= 0; vec = _mm_slli_epi64(vec, 1))
				*(uint16_t*)&OUT(rr, cc + i) = _mm_movemask_epi8(vec);
		}
	}
	if (rr == nrows) return;

	// The remainder is a block of 8x(16n+8) bits (n may be 0).
	//  Do a PAIR of 8x8 blocks in each step:
	for (cc = 0; cc <= ncols - 16; cc += 16) {
		vec = _mm_set_epi16(
			*(uint16_t const*)&INP(rr + 7, cc), *(uint16_t const*)&INP(rr + 6, cc),
			*(uint16_t const*)&INP(rr + 5, cc), *(uint16_t const*)&INP(rr + 4, cc),
			*(uint16_t const*)&INP(rr + 3, cc), *(uint16_t const*)&INP(rr + 2, cc),
			*(uint16_t const*)&INP(rr + 1, cc), *(uint16_t const*)&INP(rr + 0, cc));
		for (i = 8; --i >= 0; vec = _mm_slli_epi64(vec, 1)) {
			OUT(rr, cc + i) = h = _mm_movemask_epi8(vec);
			OUT(rr, cc + i + 8) = h >> 8;
		}
	}
	if (cc == ncols) return;

	//  Do the remaining 8x8 block:
	for (i = 0; i < 8; ++i)
		tmp.b[i] = INP(rr + i, cc);
	for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
		OUT(rr, cc + i) = _mm_movemask_epi8(tmp.x);
#undef INP
#undef OUT
}

static __m128i mm_bitshift_right(__m128i x, unsigned count)
{

	
	__m128i carry = _mm_slli_si128(x, 8);   // old compilers only have the confusingly named _mm_slli_si128 synonym
	if (count >= 64)
		return _mm_slli_epi64(carry, count - 64);  // the non-carry part is all zero, so return early
												   // else
	return _mm_or_si128(_mm_slli_epi64(x, count), _mm_srli_epi64(carry, 64 - count));



}



static __m128i mm_bitshift_left(__m128i x, unsigned count)
{

	__m128i carry = _mm_srli_si128(x, 8);   // old compilers only have the confusingly named _mm_slli_si128 synonym
	if (count >= 64)
		return _mm_srli_epi64(carry, count - 64);  // the non-carry part is all zero, so return early
	
	return _mm_or_si128(_mm_srli_epi64(x, count), _mm_slli_epi64(carry, 64 - count));
}

static __m128i mm_bitshift_right1(__m128i x, unsigned count)
{

	__m128i carry = _mm_slli_si128(x, 8);   // old compilers only have the confusingly named _mm_slli_si128 synonym
	if (count >= 64)
		return _mm_slli_epi64(carry, count - 64);  // the non-carry part is all zero, so return early
												   // else
	carry = _mm_srli_epi64(carry, 64 - count);  // After bslli shifted left by 64b

	x = _mm_slli_epi64(x, count);
	return _mm_or_si128(x, carry);
}
static __m128i mm_bitshift_left1(__m128i x, unsigned count)
{


	__m128i carry = _mm_srli_si128(x, 8);   // old compilers only have the confusingly named _mm_slli_si128 synonym
	if (count >= 64)
		return _mm_srli_epi64(carry, count - 64);  // the non-carry part is all zero, so return early
												   // else
	carry = _mm_slli_epi64(carry, 64 - count);  // After bslli shifted left by 64b

	x = _mm_srli_epi64(x, count);
	return _mm_or_si128(x, carry);
}

static void blks_bitshift_right(block* inp, int size, unsigned count)
{
	block carryOut;
	int left_shift = count % 128;
	int blk_shift = count / 128;

	block* temp = new block[size];
	memset(temp, 0, size * sizeof(block));

	//memset(out, 0, size * sizeof(block));

	if (!left_shift) //no need to shift
	{
		memcpy(inp, inp + blk_shift, (size - blk_shift) * sizeof(block));
		//memcpy(inp, temp, size * sizeof(block));
		memset(inp + (size - blk_shift), 0, blk_shift * sizeof(block));

	}
	else
	{

		int right_shift = 128 - left_shift;
		memcpy(inp, inp + blk_shift, (size - blk_shift) * sizeof(block));
		for (u64 i = 0; i < size - blk_shift - 1; i++)
		{
			inp[i] = mm_bitshift_right(inp[i], left_shift) ^ mm_bitshift_left(inp[i + 1], right_shift);
		}
		inp[size - blk_shift - 1] = mm_bitshift_right(inp[size - blk_shift - 1], left_shift);
		memset(inp + (size - blk_shift), 0, blk_shift * sizeof(block));

	}

}

static void blks_bitshift_left(block* inp, int size, unsigned count)
{
	block carryOut;
	int left_shift = count % 128;
	int blk_shift = count / 128;

	if (!left_shift) //no need to shift
	{
		memcpy(inp, inp + blk_shift, (size - blk_shift) * sizeof(block));
		memset(inp + (size - blk_shift), 0, blk_shift * sizeof(block));

	}
	else
	{

		int right_shift = 128 - left_shift;

		memcpy(inp, inp + blk_shift, (size - blk_shift) * sizeof(block));

		for (u64 i = 0; i < size - blk_shift - 1; i++)
		{
			//memcpy((u8*)&inp[i], (u8*)&(mm_bitshift_left(inp[i], left_shift) ^ mm_bitshift_right(inp[i + 1], right_shift)), sizeof(block));
			inp[i]= (mm_bitshift_left(inp[i], left_shift) ^ mm_bitshift_right(inp[i + 1], right_shift));
		}
		//memcpy(&inp[size - blk_shift - 1], &mm_bitshift_left(inp[size - blk_shift - 1], left_shift), sizeof(block));
		inp[size - blk_shift - 1]=mm_bitshift_left(inp[size - blk_shift - 1], left_shift);

	}

}

static void blks_bitshift_left(block* inp, block* out, int size, unsigned count)
{

	block carryOut;
	int left_shift = count % 128;
	int blk_shift = count / 128;

	if (!left_shift) //no need to shift
	{
		//memcpy(inp+ blk_shift, inp, (size - blk_shift) * sizeof(block));
		//memset(inp, 0, blk_shift * sizeof(block));

		memcpy(out, inp + blk_shift, (size - blk_shift) * sizeof(block));
		memset(out + (size - blk_shift), 0, blk_shift * sizeof(block));

	}
	else
	{

		int right_shift = 128 - left_shift;

		memcpy(out, inp + blk_shift, (size - blk_shift) * sizeof(block));

		for (u64 i = 0; i < size - blk_shift - 1; i++)
		{
			//	memcpy(&out[i], &(mm_bitshift_left(inp[i], left_shift) ^ mm_bitshift_right(inp[i + 1], right_shift)), sizeof(block));
			out[i]=(mm_bitshift_left(inp[i], left_shift) ^ mm_bitshift_right(inp[i + 1], right_shift));
		}
		//memcpy(&out[size - blk_shift - 1], &mm_bitshift_left(inp[size - blk_shift - 1], left_shift), sizeof(block));
		out[size - blk_shift - 1]=mm_bitshift_left(inp[size - blk_shift - 1], left_shift);

	}

}
