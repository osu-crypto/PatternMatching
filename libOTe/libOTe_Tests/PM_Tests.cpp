#include "OT_Tests.h"

#include "libOTe/TwoChooseOne/OTExtInterface.h"

#include "libOTe/Tools/Tools.h"
#include "libOTe/Tools/LinearCode.h"
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Common/Log.h>

#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"

#include "libOTe/TwoChooseOne/PM_IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/PM_IknpOtExtSender.h"

#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtReceiver.h"
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtSender.h"

#include "Common.h"
#include <thread>
#include <vector>
#include <emmintrin.h>
#include <type_traits>
#ifdef GetMessage
#undef GetMessage
#endif

#ifdef  _MSC_VER
#pragma warning(disable: 4800)
#endif //  _MSC_VER

using namespace std;
using namespace osuCrypto;



//#define SHL128(v, n) \
//({ \
//    __m128i v1, v2; \
// \
//    if ((n) >= 64) \
//    { \
//        v1 = _mm_slli_si128(v, 8); \
//        v1 = _mm_slli_epi64(v1, (n) - 64); \
//    } \
//    else \
//    { \
//        v1 = _mm_slli_epi64(v, n); \
//        v2 = _mm_slli_si128(v, 8); \
//        v2 = _mm_srli_epi64(v2, 64 - (n)); \
//        v1 = _mm_or_si128(v1, v2); \
//    } \
//    v1; \
//})

#define SHR128(v, n) \
({ \
    __m128i v1, v2; \
 \
    if ((n) >= 64) \
    { \
        v1 = _mm_srli_si128(v, 8); \
        v1 = _mm_srli_epi64(v1, (n) - 64); \
    } \
    else \
    { \
        v1 = _mm_srli_epi64(v, n); \
        v2 = _mm_srli_si128(v, 8); \
        v2 = _mm_slli_epi64(v2, 64 - (n)); \
        v1 = _mm_or_si128(v1, v2); \
    } \
    v1; \
})


namespace tests_libOTe
{
	//Modified from
	//https://mischasan.wordpress.com/2011/10/03/the-full-sse2-bit-matrix-transpose-routine/
	// with inner most loops changed to _mm_set_epi8 and _mm_set_epi16
	inline void sse_trans(uint8_t *out, uint8_t const *inp, int nrows, int ncols) {
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

	inline void sse_trans( uint8_t *inp, int nrows, int ncols) {
#   define INP(x,y) inp[(x)*ncols/8 + (y)/8]
#   define OUT(x,y) inp[(y)*nrows/8 + (x)/8]
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
	
	/*template <typename T = uint64_t>
	std::string u8_to_string(const u8 var);*/

	template <typename T = uint64_t>
	std::string u8_to_string(const u8 var) {
		std::stringstream sstr;
		const T* values = (const T*)&var;
		for (unsigned int i = 0; i < sizeof(u8) / sizeof(T); i++) {
			sstr << "0x" << std::hex << values[i] << " ";
		}
		return sstr.str();
	}

	template <typename T= uint64_t>
	std::string m128i_to_string(const __m128i var) {
		std::stringstream sstr;
		const T* values = (const T*)&var;
		for (unsigned int i = 0; i < sizeof(__m128i) / sizeof(T); i++) {
			sstr << "0x" << std::hex << values[i] << " ";
		}
		return sstr.str();
	}






	void MP_Matrix_Test_Impl() {
		

		//const uint64_t d = 1024+128;
		uint64_t length = 128 * 65536 /128;
		cout << length << endl;
		block* tt = new block[length];
		memset(tt, 0, sizeof(block) * length);
		tt[0] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
		tt[1] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
		tt[2] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
		tt[3] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
		tt[4] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
		tt[5] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
		tt[6] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
		tt[7] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);

		for (int i = 0; i <length; ++i) {
			//std::cout << m128i_to_string(tt[i]) << std::endl;
		}

		std::cout << std::endl << std::endl;
		block* tt2 = new block[128];

		sse_trans((uint8_t*)tt2, (uint8_t*)tt, 128,128);
		for (int i = 0; i <128; ++i) {
			//std::cout << tt2[i] << std::endl;
		std::cout << m128i_to_string(tt2[i]) << std::endl;
		}

		/*for (auto& d : tt)
		{
			if (neq(d, _mm_set_epi64x(0, 0xFF)))
			{
				std::cout << "expected" << std::endl;
				std::cout << _mm_set_epi64x(0xF, 0) << std::endl << std::endl;

				printMtx(tt);

				throw UnitTestFail();
			}
		}*/

	}


   
    void printMtx(std::array<block, 128>& data)
    {
        for (auto& d : data)
        {
            std::cout << d << std::endl;
        }
    }

    void PM_Transpose_Test_Impl()
    {
        {

            std::array<block, 128> data;
            memset((u8*)data.data(), 0, sizeof(data));

            data[0] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[1] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[2] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[3] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[4] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[5] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[6] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[7] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);

            //printMtx(data);
            eklundh_transpose128(data);



            for (auto& d : data)
            {
                if (neq(d, _mm_set_epi64x(0, 0xFF)))
                {
                    std::cout << "expected" << std::endl;
                    std::cout << _mm_set_epi64x(0xF, 0) << std::endl << std::endl;

                    printMtx(data);

                    throw UnitTestFail();
                }
            }
        }
        {


            std::array<block, 128> data;
            memset((u8*)data.data(), 0, sizeof(data));

            data[0] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[1] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[2] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[3] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[4] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[5] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[6] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
            data[7] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);

            sse_transpose128(data);


            for (auto& d : data)
            {
                if (neq(d, _mm_set_epi64x(0, 0xFF)))
                {
                    std::cout << "expected" << std::endl;
                    std::cout << _mm_set_epi64x(0xF, 0) << std::endl << std::endl;

                    printMtx(data);

                    throw UnitTestFail();
                }
            }
        }

        {
            PRNG prng(ZeroBlock);

            std::array<std::array<block, 8>, 128> data;

            prng.get((u8*)data.data(), sizeof(block) * 8 * 128);


            std::array<std::array<block, 8>, 128> data2 = data;

            sse_transpose128x1024(data);


            for (u64 i = 0; i < 8; ++i)
            {

                std::array<block, 128> sub;

                for (u64 j = 0; j < 128; ++j)
                {
                    sub[j] = data2[j][i];
                }

                sse_transpose128(sub);

                for (u64 j = 0; j < 128; ++j)
                {
                    if (neq(sub[j], data[j][i]))
                        throw UnitTestFail();
                }
            }

        }
    }

	void OT_100Receive_Test(BitVector& choiceBits, gsl::span<block> recv, gsl::span<std::array<block, 2>>  sender)
	{

		for (u64 i = 0; i < choiceBits.size(); ++i)
		{

			u8 choice = choiceBits[i];
			const block & revcBlock = recv[i];
			//(i, choice, revcBlock);
			const block& senderBlock = sender[i][choice];

			//if (i%512==0) {
			//    std::cout << "[" << i << ",0]--" << sender[i][0] << std::endl;
			//    std::cout << "[" << i << ",1]--" << sender[i][1] << std::endl;
			//    std::cout << (int)choice << "-- " << recv[i] << std::endl;
			//}
			if (neq(revcBlock, senderBlock))
				throw UnitTestFail();

			if (eq(revcBlock, sender[i][1 ^ choice]))
				throw UnitTestFail();
		}

	}

    void PM_TransposeMatrixView_Test_Impl()
    {
		
        {

            PRNG prng(ZeroBlock);

            std::array<block, 128> data;
            prng.get(data.data(), data.size());
            std::array<block, 128> data2;

            MatrixView<block> dataView(data.begin(), data.end(), 1);
            MatrixView<block> data2View(data2.begin(), data2.end(), 1);

            sse_transpose(dataView, data2View);

            sse_transpose128(data);




            for (u64 i = 0; i < 128; ++i)
            {
                if (neq(data[i], data2[i]))
                {
                    std::cout << i << "\n";
                    printMtx(data);
                    std::cout << "\n";
                    printMtx(data2);

                    throw UnitTestFail();
                }
            }
        }


        {
            PRNG prng(ZeroBlock);

            std::array<std::array<block, 8>, 128> data;

            prng.get((u8*)data.data(), sizeof(block) * 8 * 128);


            std::array<std::array<block, 8>, 128> data2;

            MatrixView<block> dataView((block*)data.data(), 128, 8);
            MatrixView<block> data2View((block*)data2.data(), 128 * 8, 1);
            sse_transpose(dataView, data2View);


            for (u64 i = 0; i < 8; ++i)
            {
                std::array<block, 128> data128;

                for (u64 j = 0; j < 128; ++j)
                {
                    data128[j] = data[j][i];
                }

                sse_transpose128(data128);


                for (u64 j = 0; j < 128; ++j)
                {
                    if (neq(data128[j], data2View[i * 128 + j][0]))
                        throw UnitTestFail();
                }
            }

        }


        {
            PRNG prng(ZeroBlock);

            //std::array<std::array<std::array<block, 8>, 128>, 2> data;

            Matrix<block> dataView(208, 8);
            prng.get((u8*)dataView.data(), sizeof(block) *dataView.bounds()[0] * dataView.stride());

            Matrix<block> data2View(1024, 2);
            memset(data2View.data(), 0, data2View.bounds()[0] * data2View.stride() * sizeof(block));
            sse_transpose(dataView, data2View);

            for (u64 b = 0; b < 2; ++b)
            {

                for (u64 i = 0; i < 8; ++i)
                {
                    std::array<block, 128> data128;

                    for (u64 j = 0; j < 128; ++j)
                    {
                        if (dataView.bounds()[0] > 128 * b + j)
                            data128[j] = dataView[128 * b + j][i];
                        else
                            data128[j] = ZeroBlock;
                    }

                    sse_transpose128(data128);

                    for (u64 j = 0; j < 128; ++j)
                    {
                        if (neq(data128[j], data2View[i * 128 + j][b]))
                        {
                            std::cout << "failed " << i << "  " << j << "  " << b << std::endl;
                            std::cout << "exp: " << data128[j] << "\nact: " << data2View[i * 128 + j][b] << std::endl;
                            throw UnitTestFail();
                        }
                    }
                }
            }
        }

        {
            PRNG prng(ZeroBlock);

            Matrix<u8> in(16, 8);
            prng.get((u8*)in.data(), sizeof(u8) *in.bounds()[0] * in.stride());

            Matrix<u8> out(63, 2);
            sse_transpose(in, out);


            Matrix<u8> out2(64, 2);
            sse_transpose(in, out2);

            for (u64 i = 0; i < out.bounds()[0]; ++i)
            {
                if (memcmp(out[i].data(), out2[i].data(), out[i].size()))
                {
                    std::cout << "bad " << i << std::endl;
                    throw UnitTestFail();
                }
            }
        }

        {
            PRNG prng(ZeroBlock);

            //std::array<std::array<std::array<block, 8>, 128>, 2> data;

            Matrix<u8> in(25, 9);
            Matrix<u8> in2(32, 9);

            prng.get((u8*)in.data(), sizeof(u8) *in.bounds()[0] * in.stride());
            memset(in2.data(), 0, in2.bounds()[0] * in2.stride());

            for (u64 i = 0; i < in.bounds()[0]; ++i)
            {
                for (u64 j = 0; j < in.stride(); ++j)
                {
                    in2[i][j] = in[i][j];
                }
            }

            Matrix<u8> out(72, 4);
            Matrix<u8> out2(72, 4);

            sse_transpose(in, out);
            sse_transpose(in2, out2);

            for (u64 i = 0; i < out.bounds()[0]; ++i)
            {
                for (u64 j = 0; j < out.stride(); ++j)
                {
                    if (out[i][j] != out2[i][j])
                    {
                        std::cout << (u32)out[i][j] << " != " << (u32)out2[i][j] << std::endl;
                        throw UnitTestFail();
                    }
                }
            }
        }
    }


	//template<typename T>
	//typename std::enable_if_t<std::is_pod<T>::value, void>
	//	Channel::asyncSend(const T * buffT, u64 sizeT)
	//{
	//	u8* buff = (u8*)buffT;
	//	auto size = sizeT * sizeof(T);

	//	// not zero and less that 32 bits
	//	Expects(size - 1 < u32(-2) && mBase->mSendStatus == Status::Normal);

	//	auto op = std::unique_ptr<IOOperation>(new PointerSizeBuff(buff, size, IOOperation::Type::SendData));
	//	mBase->getIOService().dispatch(mBase.get(), std::move(op));
	//}


	__m128i SHL128(__m128i v, int n)
	{
		__m128i v1, v2;

		if ((n) >= 64)
		{
			v1 = _mm_slli_si128(v, 8);
			v1 = _mm_slli_epi64(v1, (n)-64);
		}
		else
		{
			v1 = _mm_slli_epi64(v, n);
			v2 = _mm_slli_si128(v, 8);
			v2 = _mm_srli_epi64(v2, 64 - (n));
			v1 = _mm_or_si128(v1, v2);
		}
		return  v1;
	}

	//void shift()


	// doesn't optimize for the special-case where count%8 = 0
	// could maybe do that in gcc with if(__builtin_constant_p(count)) { if (!count%8) return ...; }



	//----------------------------------------------------------------------------
	// bit shift left a 128-bit value using xmm registers
	//          __m128i *data - data to shift
	//          int count     - number of bits to shift
	// return:  __m128i       - carry out bit(s)

	static __m128i bitShiftLeft128xmm(__m128i *data, int count)
	{
		__m128i innerCarry, carryOut;

		innerCarry = _mm_srli_epi64(*data, 64 - count);      // carry outs in bit 0 of each qword
		carryOut = _mm_shuffle_epi32(innerCarry, 0xFE);    // upper carry in xmm bit 0, others zero
		innerCarry = _mm_shuffle_epi32(innerCarry, 0xCF);    // lower carry in xmm bit 64, others zero
		*data = _mm_slli_epi64(*data, count);                // shift all qwords left
		*data = _mm_or_si128(*data, innerCarry);             // propagate carry out from low qword
		return carryOut;
	}

	//----------------------------------------------------------------------------
	// bit shift left a 256-bit value using xmm registers
	//          __m128i *data - data to shift, ls part stored first 
	//          int count     - number of bits to shift
	// return:  __m128i       - carry out bit(s)

	static __m128i bitShiftLeft256xmm(__m128i *data, int count)
	{
		__m128i carryOut0, carryOut1;

		carryOut0 = bitShiftLeft128xmm(&data[0], count);
		carryOut1 = bitShiftLeft128xmm(&data[1], count);
		data[1] = _mm_or_si128(data[1], carryOut0);
		return carryOut1;
	}

	
	//----------------------------------------------------------------------------
	//
	// shiftLeftOnce - shift left extended integer once using SSE registers
	//
	static void shiftLeftOnceSse(block *source, block *dest, int activeBits)
	{
		int sourceIndex, uint128Count;
		__m128i current, next, shiftedBits, recoveredBits, carryOut;

		uint128Count = activeBits / 128;
		sourceIndex = uint128Count - 1;
		next = source[sourceIndex];
		while (sourceIndex)
		{
			current = source[sourceIndex];
			next = source[sourceIndex - 1];
			carryOut = _mm_srli_epi64(next, 63);                    // bit 64 is bit 127 from next part
			carryOut = _mm_shuffle_epi32(carryOut, 0xFE);           // bit 0 is bit 127 from next part, other bits are zero
			recoveredBits = _mm_srli_epi64(current, 63);            // bit 0 is bit 63 of current part 
			recoveredBits = _mm_shuffle_epi32(recoveredBits, 0xCF); // bit 64 is bit 63 of current part, other bits are zero
			recoveredBits = _mm_or_si128(recoveredBits, carryOut);  // bits 0 and 64 are bits lost by the SSE shift
			shiftedBits = _mm_slli_epi64(current, 1);
			dest[sourceIndex] = _mm_or_si128(shiftedBits, recoveredBits);
			current = next;
			sourceIndex--;
		}

		// the final shift zero fills
		dest[0] = _mm_or_si128(_mm_slli_epi64(next, 1), _mm_shuffle_epi32(_mm_srli_epi64(next, 63), 0xCF));
	}

	template <typename T>
	constexpr void shift_array_left(T *arr, const size_t size, const size_t bits, const bool zero = false) {
		const size_t chunks = bits / (8 * sizeof(T));

		if (chunks >= size) {
			if (zero) {
				memset(arr, 0, size);
			}
			return;
		}

		if (chunks) {
			memmove(arr, arr + chunks, size - chunks);
			if (zero) {
				memset(arr + size - chunks, 0, chunks);
			}
		}

		const size_t left = bits % (8 * sizeof(T));

		// If we have non directly addressable bits left we need to move the whole thing one by one.
		if (left) {
			const size_t right = (8 * sizeof(T)) - left;
			const size_t l = size - chunks - 1;
			for (size_t i = 0; i < l; i++) {
				arr[i] = ((arr[i] << left) | (arr[i + 1] >> right));
			}
			arr[l] = (arr[l] << left) ;
		}
	}

	/*template <typename T = u8>
	constexpr void shift_array_left(T *arr, const size_t size, const size_t bits, const bool zero = false)*/;

	void blks_shift_left(block* &arr, const size_t size, const size_t bits, const bool zero = false) {
	
		u8* tmp = new u8[sizeof(block) * size];
		memcpy(tmp, arr, sizeof(block) * size);

		shift_array_left(tmp, sizeof(block) * size, bits,zero);

		for (u16 i = 0; i < size; i++)
			memcpy(&arr[i], tmp+i* sizeof(block), sizeof(block));
	}
	

	static __m128i mm_bitshift_left(__m128i x, unsigned count)
	{
		__m128i carry = _mm_slli_si128(x, 8);   // old compilers only have the confusingly named _mm_slli_si128 synonym
		if (count >= 64)
			return _mm_slli_epi64(carry, count - 64);  // the non-carry part is all zero, so return early
													   // else
		carry = _mm_srli_epi64(carry, 64 - count);  // After bslli shifted left by 64b

		x = _mm_slli_epi64(x, count);
		return _mm_or_si128(x, carry);
	}

	static __m128i mm_bitshift_right(__m128i x, unsigned count)
	{
		__m128i carry = _mm_srli_si128(x, 8);   // old compilers only have the confusingly named _mm_slli_si128 synonym
		if (count >= 64)
			return _mm_srli_epi64(carry, count - 64);  // the non-carry part is all zero, so return early
													   // else
		carry = _mm_slli_epi64(carry, 64 - count);  // After bslli shifted left by 64b

		x = _mm_srli_epi64(x, count);
		return _mm_or_si128(x, carry);
	}

	static void blks_bitshift_left(block* inp, int size, unsigned count)
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
				inp[i] = mm_bitshift_left(inp[i], left_shift) ^ mm_bitshift_right(inp[i + 1], right_shift);
			}
			inp[size - blk_shift - 1] = mm_bitshift_left(inp[size - blk_shift - 1], left_shift);
			memset(inp + (size - blk_shift), 0, blk_shift * sizeof(block));

		}

	}

	/*static void blks_bitshift_left(block*& blk, int size, unsigned count)
	{
		block* tmp = new block[size];
		memcpy(blk,tmp,size*)
		blks_bitshift_left(blk, blk, size, count);
	}*/
	void Shift_Test_Impl()
	{

		__m128i va = _mm_setr_epi8(0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f);
	

		u64 numBlk = 4;
		block* sendMsg1 = new block[numBlk];
		block* out = new block[numBlk];
		sendMsg1[0] = _mm_setr_epi8(0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f);
		sendMsg1[1] = _mm_setr_epi8(0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f);

		


		sendMsg1[2] = _mm_set_epi64x(1, 0xFFFFFFFFFFFFF000);
		sendMsg1[3] =_mm_set_epi64x(0xFFFFFFFFFFFFF000, 1);
		std::cout << sendMsg1[0] << " ";
		std::cout << sendMsg1[1] << endl;

		for (size_t i = 0; i < 128*3; i+=4)
		{
			blks_bitshift_left(sendMsg1,  4, 4);
			std::cout << sendMsg1[0] << " ";
			std::cout << sendMsg1[1] << " ";
			std::cout << sendMsg1[2] << " ";
			std::cout << sendMsg1[3] << endl;
		}

		//block* a = new block[numBlk];

		block a;
		for (size_t i = 0; i < 15; i++)
		{
			//block a;
			a = mm_bitshift_left(sendMsg1[1], i);
			//std::cout << a << endl;
			//a=mm_bitshift_right(sendMsg1[0],i);
			//std::cout << a << endl;

		}

		//std::cout << endl;
		//std::cout  << endl;

		for (size_t i = 0; i < 15; i++)
		{
			//block a;
			a = mm_bitshift_right(sendMsg1[0], i);
			a = mm_bitshift_left(a, i);
			//std::cout << a << endl;
			//a=mm_bitshift_right(sendMsg1[0],i);
			//std::cout << a << endl;

		}
		


		//int bits = 8;
		//const size_t left = bits % 8;

		//// If we have non directly addressable bits left we need to move the whole thing one by one.
		//if (left) {
		///*	const size_t right = (8 * sizeof(T)) - left;
		//	const size_t l = size - chunks - 1;
		//	for (size_t i = 0; i < l; i++) {
		//		arr[i] = ((arr[i] << left) & ~left) | (arr[i + 1] >> right);
		//	}
		//	arr[l] = (arr[l] << left) & ~left;*/
	
		//	for (size_t i = 0; i < left; i++)
		//	{

		//	}
		//}
		//else {
		//	const size_t chunks = bits / 8;
		//	memmove((u8*)&sendMsg1, (u8*)&sendMsg1 + chunks, sizeof(block)*numBlk - chunks);
		//}

		//

		//bitShiftLeft256xmm(sendMsg1, 4);
		//blks_shift_left(sendMsg1, 2, 16,true);

		
		//std::cout << sendMsg1[0] << " ";
		//std::cout << sendMsg1[1] << endl;
	}


	void PM_IknpOtExt_100Receive_Test_Impl()
	{
#pragma region channel
		setThreadName("Sender");
		IOService ios(0);
		Endpoint ep0(ios, "127.0.0.1", 1212, EpMode::Server, "ep");
		Endpoint ep1(ios, "127.0.0.1", 1212, EpMode::Client, "ep");
		Channel senderChannel = ep1.addChannel("chl", "chl");
		Channel recvChannel = ep0.addChannel("chl", "chl");
#pragma endregion

		PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
		PRNG prng1(_mm_set_epi32(4253233465, 334565, 0, 235));

		u64 numBlkText = 5;
		u64 numBlkPattern = 2;

		u64 numOTs = numBlkPattern*128;

		std::vector<block> recvSeeds(numOTs), baseRecv(128);
		std::vector<std::array<block, 2>> sendSeeds(numOTs), baseSend(128);
		BitVector choices(numOTs), baseChoice(128);
		choices.randomize(prng0);
		baseChoice.randomize(prng0);
		
		for (size_t i = 0; i < numOTs; i++)
		{
			choices[i] = 1;
		}

		int blk_step = 1;
		block* pmRecvMsg = new block[numOTs*numBlkText];
		block* pmSendMsg =new block[numOTs*numBlkText];

		block* trans_pmRecvMsg = new block[numOTs*blk_step];
		block* trans_pmSendMsg = new block[numOTs*blk_step];


		
		block* delta_off = new block[numBlkText];
		block* text = new block[numBlkText];
		prng0.get(delta_off, numBlkText);
		prng0.get(text, numBlkText);

		block* pattern = new block[numBlkText];
//		block* text = new block[numBlkText];
		prng0.get(delta_off, numBlkText);
	//	prng0.get(text, numBlkText);



		for (size_t i = 0; i < numBlkText; i++)
		{
			//delta_off[i] = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
			//delta_off[i] = _mm_set_epi64x(0,0);
			//delta_off[i] = AllOneBlock;
		}
		cout << delta_off[0] << endl;

		
		int rowId = 129;
		int colId = 2;
		int testId = rowId*numBlkText+colId;

		std::cout << IoStream::lock;
		cout << "choice[j] " << choices[rowId] << endl;
		std::cout << IoStream::unlock;
		//cout << delta_off[0] << endl;
		
		block* tmp_s = new block[numOTs*blk_step];
		block* tmp_r = new block[numOTs*blk_step];

		//for (u64 i = 0; i < numOTs; i++)
		//{
		//	pmRecvMsg[i] = new block[numBlkText];
		//	pmSendMsg[i] = new block[numBlkText];
		//}
		

		/*for (u64 i = 0; i < 128; ++i)
		{
			baseRecv[i] = baseSend[i][baseChoice[i]];
		}*/

		PM_IknpOtExtSender sender;
		PM_IknpOtExtReceiver recv;

		std::thread thrd = std::thread([&]() {
			
			recv.setBaseOts(baseSend);
			recv.receive(choices, recvSeeds, prng0, recvChannel);
			std::cout << recvSeeds[rowId] << std::endl;


			//matrix |p|x|t|
			
			//bool = true;


			for (u64 i = 0; i < numOTs; ++i)
			{
				PRNG prg_pm;
				prg_pm.SetSeed(recvSeeds[i]);

				//block* test_r = new block[numBlkText];
				//prg_pm.get(test_r, numBlkText);
				prg_pm.get(pmRecvMsg+i*numBlkText, numBlkText);

				block* tmp_recv = new block[numBlkText];
				recvChannel.recv((u8*)tmp_recv, sizeof(block)*(numBlkText));
				//memset((u8*)tmp_recv, 0, numBlkText * sizeof(block));

				if (i == rowId) {				
				std::cout << IoStream::lock;
				cout << "r tmp_recv[" << i << "] " << tmp_recv[0] << endl;
				std::cout << IoStream::unlock;

				}

				if (choices[i] == 1)
					for (u64 j = 0; j < numBlkText; ++j)					
					{
						pmRecvMsg[i*numBlkText + j] = pmRecvMsg[i*numBlkText + j] ^ tmp_recv[j];
						//cout << "pmRecvMsg[" << i << "]["<<j<<"] " << pmRecvMsg[i*numBlkText + j] << endl;
					}
					/*else
					{
						pmRecvMsg[i*numBlkText + j] = test_r[j];
					}*/
			}

			std::cout << IoStream::lock;
			cout << "pmRecvMsg[j] " << pmRecvMsg[testId] << endl;
			std::cout << IoStream::unlock;


			//shift
			for (u64 i = 1; i < numOTs; ++i)
				blks_bitshift_left(pmRecvMsg + i*numBlkText, numBlkText, i);

			std::cout << IoStream::lock;
			cout << "shift pmRecvMsg[j] " << pmRecvMsg[testId] << endl;
			std::cout << IoStream::unlock;


			
			for (size_t i = 0; i < numOTs; i++)
			{
				memcpy(tmp_r + i*blk_step, pmRecvMsg + i*numBlkText, blk_step * sizeof(block));
				std::cout << IoStream::lock;
				//cout << "tmp_r["<<i<<"] " << tmp_r[i] << endl;
				std::cout << IoStream::unlock;

			}
			
#if 1
			sse_trans((uint8_t*)trans_pmRecvMsg, (uint8_t*)tmp_r, numOTs, 128*blk_step);

			std::cout << IoStream::lock;
			cout << "trans_pmRecvMsg[j] " << trans_pmRecvMsg[0] << endl;
			std::cout << IoStream::unlock;
#endif

		});


		sender.setBaseOts(baseRecv, baseChoice);
		sender.send(sendSeeds, prng1, senderChannel);

		std::cout << sendSeeds[rowId][0] << std::endl;
		std::cout << sendSeeds[rowId][1] << std::endl;


		//matrix |p|x|t|
		for (u64 i = 0; i < numOTs; ++i)
		{
			std::vector<PRNG> prg_pm(2);
			block* delta_send = new block[numBlkText];

			prg_pm[0].SetSeed(sendSeeds[i][0]);
			prg_pm[1].SetSeed(sendSeeds[i][1]);

			//std::vector<block*> test_s(2);			
			//test_s[0] = new block[numBlkText];
			//test_s[1] = new block[numBlkText];

			//prg_pm[0].get(test_s[0], numBlkText);
			//prg_pm[1].get(test_s[1], numBlkText);

			for (u64 j = 0; j < numBlkText; ++j)
			{				
				/*pmSendMsg[i*numBlkText + j] = test_s[0][j] ^ delta_off[j];
				delta_send[j] = pmSendMsg[i*numBlkText + j] ^ test_s[1][j];*/
				pmSendMsg[i*numBlkText + j] = prg_pm[0].get<block>() ^ delta_off[j];
				delta_send[j] = pmSendMsg[i*numBlkText + j] ^ prg_pm[1].get<block>();
			}
			if (i == rowId)
			{
				std::cout << IoStream::lock;
			cout << "s tmp_s_recv["<<i<<"] " << delta_send[0] << endl;
			std::cout << IoStream::unlock;
			}
			senderChannel.asyncSend((u8*)delta_send, sizeof(block)*(numBlkText));
		}

		std::cout << IoStream::lock;
		cout << "pmSendMsg[j] " << pmSendMsg[testId] << endl;
		std::cout << IoStream::unlock;

		
		//shift
		for (u64 i = 1; i < numOTs; ++i)
			blks_bitshift_left(pmSendMsg + i*numBlkText, numBlkText, i);

		std::cout << IoStream::lock;
		cout << "shift pmSendMsg[j] " << pmSendMsg[testId] << endl;
		std::cout << IoStream::unlock;

		//memset((uint8_t*)pmSendMsg, 0, sizeof(block)*numBlkText*numOTs);

		
		for (size_t i = 0; i < numOTs; i++)
		{
			//for (size_t j = 0; j < blk_step; j++)
				//tmp[i*blk_step + j] = pmSendMsg[i*numBlkText + j];
				memcpy(tmp_s +i*blk_step,pmSendMsg+i*numBlkText, blk_step*sizeof(block));
				std::cout << IoStream::lock;
				//cout << "tmp_s["<<i<<"] " << tmp_s[i] << endl;
				std::cout << IoStream::unlock;
		}

#if 1

		sse_trans((uint8_t*)trans_pmSendMsg, (uint8_t*)tmp_s, numOTs, 128*blk_step);
	
		std::cout << IoStream::lock;
		cout << "trans_pmSendMsg[j] " << trans_pmSendMsg[0] << endl;
		std::cout << IoStream::unlock;
#endif

		thrd.join();

		//check seeds
		for (u64 i = 0; i < choices.size(); ++i)
		{
			u8 choice = choices[i];
			const block & revcBlock = recvSeeds[i];
			const block& senderBlock = sendSeeds[i][choice];

			if (neq(revcBlock, senderBlock))
						cout << i << ": " << choices[i] << " " <<senderBlock << "  " << revcBlock  << endl;

			if (eq(revcBlock, sendSeeds[i][1 ^ choice]))
				cout << i << ": " << choices[i] << " " << senderBlock << "  " << sendSeeds[i][1 ^ choice] << endl;
		}

		//check OT
		for (size_t i = 0; i < numOTs*blk_step; i++)
		{
			const block & revcBlock = pmRecvMsg[i*numBlkText];
			const block& senderBlock = pmSendMsg[i*numBlkText];
			if(choices[i]==1)
			if (neq(revcBlock, senderBlock))
			{	//throw UnitTestFail();
				//cout << i << ": " << tmp_r[i] << "  " << tmp_s[i] << "\t ";
				block a = senderBlock^delta_off[0];
				cout << i << ": " << choices[i] << " " << senderBlock << "  " << revcBlock << endl;
				cout << "\t " << sendSeeds[i][choices[i]] << "  " << recvSeeds[i] << endl;
				cout << "\t " << a<< endl;

				PRNG prg_r;
				std::vector<PRNG> prg_s(2);
				block* test_r = new block[numBlkText];
				std::vector<block*> test_s(2);
				test_s[0] = new block[numBlkText];
				test_s[1] = new block[numBlkText];

				prg_s[0].SetSeed(sendSeeds[i][0]);
				prg_s[1].SetSeed(sendSeeds[i][1]);
				prg_r.SetSeed(recvSeeds[i]);

				prg_r.get(test_r, numBlkText);
				prg_s[0].get(test_s[0], numBlkText);
				prg_s[1].get(test_s[1], numBlkText);
				for (size_t j = 0; j < numBlkText; j++)
				{
					if (neq(test_r[j], test_s[choices[i]][j]))
					{
						cout << "test\t " << j << ": " << test_r[j] << " " << test_s[choices[i]][j] << endl;

					}
				}

				block* dtest_r = new block[numBlkText];
				block* dtest_s = new block[numBlkText];
				//std::vector<block*> dtest_s(2);
				//dtest_s[0] = new block[numBlkText];
				//dtest_s[1] = new block[numBlkText];

				for (u64 j = 0; j < numBlkText; ++j)
				{
					//dtest_s[j] = test_s[0][j] ^ delta_off[j];
					block del_send = test_s[0][j] ^ delta_off[j] ^ test_s[1][j];
					if(choices[i]==1)
						dtest_r[j] = test_r[j] ^ del_send;
					else 
						dtest_r[j] = test_r[j] ;

				}

				for (size_t j = 0; j < numBlkText; j++)
				{
					if (choices[i] == 1)
						if (neq(test_s[0][j]^ delta_off[j], dtest_r[j]))
						{
							cout << "dtest\t " << j << ": " << dtest_s[j] << " " << dtest_s[j] << endl;
						}

					if (choices[i] == 0)
						if (neq(test_s[0][j], dtest_r[j]))
						{
							cout << "dtest\t " << j << ": " << dtest_r[j] << " " << test_s[0][j] << endl;
						}
				}

				/*if (choices[i])
					for (u64 j = 0; j < numBlkText; ++j)
					{
						pmRecvMsg[i*numBlkText + j] = pmRecvMsg[i*numBlkText + j] ^ tmp_recv[j];
					}*/

			}
		}

		//check transpose

		senderChannel.close();
		recvChannel.close();


		ep1.stop();
		ep0.stop();

		ios.stop();

		//senderNetMgr.Stop();
		//recvNetMg
	}



    void PM_IknpOtExt_100Receive_Test_Impl1()
    {
        setThreadName("Sender");

        IOService ios(0);
        Endpoint ep0(ios, "127.0.0.1", 1212, EpMode::Server, "ep");
        Endpoint ep1(ios, "127.0.0.1", 1212, EpMode::Client, "ep");
        Channel senderChannel = ep1.addChannel("chl", "chl");
        Channel recvChannel = ep0.addChannel("chl", "chl");

        PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
        PRNG prng1(_mm_set_epi32(4253233465, 334565, 0, 235));

        u64 numOTs = 200;

        std::vector<block> recvSeeds(numOTs), baseRecv(128);
        std::vector<std::array<block, 2>> sendSeeds(numOTs), baseSend(128);
        BitVector choices(numOTs), baseChoice(128);
        choices.randomize(prng0);
        baseChoice.randomize(prng0);

		std::vector<block*> pmRecvMsg(numOTs);
		std::vector<block*> pmTempRecvMsg(numOTs);
		std::vector<std::array<block*, 2>> pmSendMsg(numOTs);
		u64 numBlkText = 4;
		u64 numBlkPattern= 1;

		

		for (u64 i = 0; i < numOTs; i++)
		{
			pmRecvMsg[i] = new block[numBlkText];
			pmTempRecvMsg[i] = new block[numBlkText];
			pmSendMsg[i][0] = new block[numBlkText];
			pmSendMsg[i][1] = new block[numBlkText];
			for (u64 j = 0; j < numBlkText; j++)
			{
				pmSendMsg[i][0][j] = _mm_set_epi64x(1, 1); //prng0.get<block>();
				pmSendMsg[i][1][j] = _mm_set_epi64x(1, 1); //prng0.get<block>();
			}

		}
		//prng0.get((u8*)pmSendMsg.data()->data(), numBlkText * sizeof(block) * 2 * pmSendMsg.size());

		std::cout << "pmSendMsg[0][0][0]: " <<pmSendMsg[0][0][0] << std::endl;
		std::cout << "pmSendMsg[0][1][0]: " << pmSendMsg[0][1][0] << std::endl;


       for (u64 i = 0; i < 128; ++i)
        {
            baseRecv[i] = baseSend[i][baseChoice[i]];
        }

        PM_IknpOtExtSender sender;
		PM_IknpOtExtReceiver recv;

        std::thread thrd = std::thread([&]() {
			recv.setBaseOts(baseSend);
            recv.receive(choices, recvSeeds, prng0, recvChannel);
			std::cout << recvSeeds[0] << std::endl;

			PRNG prg_pm;
			for (u64 i = 0; i < numOTs; ++i)
			{
				prg_pm.SetSeed(recvSeeds[i]);
				if (choices[i])
				{
					recvChannel.recv((u8*)&pmTempRecvMsg[i], sizeof(block)*(numBlkText - numBlkPattern));
					recvChannel.recv((u8*)&pmRecvMsg[i], sizeof(block)*(numBlkText - numBlkPattern));

				/*	if (i == 0)
						std::cout << "pmRecvMsg[i][0]: " << pmRecvMsg[i][0] << std::endl;*/
				}
				else
				{
					recvChannel.recv((u8*)&pmRecvMsg[i], sizeof(block)*(numBlkText - numBlkPattern));
					recvChannel.recv((u8*)&pmTempRecvMsg[i], sizeof(block)*(numBlkText - numBlkPattern));
				}

				if (i == 0)
					std::cout << "r pmRecvMsg[i][0]: " << pmRecvMsg[i][0] << std::endl;


				for (u64 j = 0; j < numBlkText - numBlkPattern; j++)
				{
					block tmp = prg_pm.get<block>();
					if (i == 0 && j == 0)
						std::cout << "tmp " << tmp << std::endl;
					pmRecvMsg[i][j] = pmRecvMsg[i][j] ^ tmp;
				}
				if (i == 0)
					std::cout << "r pmRecvMsg[i][0]: " << pmRecvMsg[i][0] << std::endl;

			}
		
			

        });



        //{
        //    std::lock_guard<std::mutex> lock(Log::mMtx);
        //    for (u64 i = 0; i < baseOTs.receiver_outputs.size(); ++i)
        //    {
        //        std::cout << "i  " << baseOTs.receiver_outputs[i] << " " << (int)baseOTs.receiver_inputs[i] << std::endl;
        //    }
        //}
        sender.setBaseOts(baseRecv, baseChoice);
        sender.send(sendSeeds, prng1, senderChannel);
	
		std::cout << sendSeeds[0][0] << std::endl;
		std::cout << sendSeeds[0][1] << std::endl;
/*
		std::cout << sendSeeds[0][0] << std::endl;
		std::cout << sendSeeds[0][1] << std::endl;*/

		PRNG prg_pm1;
		PRNG prg_pm2;
		for (u64 i = 0; i < numOTs; ++i)
		{
			prg_pm1.SetSeed(sendSeeds[i][0]);
			prg_pm2.SetSeed(sendSeeds[i][1]);


			if (i == 0)
			{
				std::cout << "bs pmSendMsg[i][0]: " << pmSendMsg[i][0][0] << std::endl;
				std::cout << "bs pmSendMsg[i][1]: " << pmSendMsg[i][1][0] << std::endl;
			}

			for (u64 j = 0; j < numBlkText; j++)
			{
				block tmp1 = prg_pm1.get<block>();
				block tmp2 = prg_pm2.get<block>();
				if (i == 0 && j == 0)
				{
					std::cout << "tmp1 " << tmp1 << std::endl;
					std::cout << "tmp2 " << tmp2 << std::endl;
				}
				pmSendMsg[i][0][j] = pmSendMsg[i][0][j] ^ tmp1;
				pmSendMsg[i][1][j] = pmSendMsg[i][1][j] ^ tmp2;
			}

			u8* C = new u8[sizeof(block) * numBlkText];
			memcpy(C, pmSendMsg[i][0], sizeof(block) * numBlkText);
			shift_array_left(C, sizeof(block) * numBlkText, i);

			if (i == 0)
			{
				block* sendMsg2 = new block[numBlkText];
				for (size_t k = 0; k < numBlkText; k++)
				{
					memcpy(&sendMsg2[k], C+k* sizeof(block), sizeof(block));
					cout << "sendMsg2[k]" << sendMsg2[k] << endl;
				}
			}



			senderChannel.asyncSend(C, sizeof(block)*(numBlkText - numBlkPattern));

			//senderChannel.asyncSend((u8*)&pmSendMsg[i][0]+1, sizeof(block)*(numBlkText- numBlkPattern));
			
			//u8* C = new u8[sizeof(block) * numBlkText];
			memcpy(C, pmSendMsg[i][1], sizeof(block) * numBlkText);
			shift_array_left(C, sizeof(block) * numBlkText, i);

			if (i == 0)
			{
				block* sendMsg2 = new block[numBlkText];
				for (size_t k = 0; k < numBlkText; k++)
				{
					memcpy(&sendMsg2[k], C + k * sizeof(block), sizeof(block));
					cout << "sendMsg2[k]" << sendMsg2[k] << endl;
				}
			}

			senderChannel.asyncSend(C, sizeof(block)*(numBlkText - numBlkPattern));


			
		}



		thrd.join();

        //for (u64 i = 0; i < baseOTs.receiver_outputs.size(); ++i)
        //{
        //    std::cout << sender.GetMessage(i, 0) << " " << sender.GetMessage(i, 1) << "\n" << recv.GetMessage(1) << "  " << recv.mChoiceBits[i] << std::endl;
        //}
        OT_100Receive_Test(choices, recvSeeds, sendSeeds);




        senderChannel.close();
        recvChannel.close();


        ep1.stop();
        ep0.stop();

        ios.stop();

        //senderNetMgr.Stop();
        //recvNetMg
    }

	


	void PM_KkrtNcoOt_Test_Impl()
	{
		setThreadName("Sender");

		PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
		PRNG prng1(_mm_set_epi32(4253465, 3434565, 234435, 23987025));

		// The total number that we wish to do
		u64 numOTs = 1030;

		KkrtNcoOtSender sender;
		KkrtNcoOtReceiver recv;

		// get up the parameters and get some information back. 
		//  1) false = semi-honest
		//  2) 40  =  statistical security param.
		//  3) numOTs = number of OTs that we will perform
		sender.configure(false, 40, 128);
		recv.configure(false, 40, 128);

		// the number of base OT that need to be done
		u64 baseCount = sender.getBaseOTCount();

		// Fake some base OTs
		std::vector<block> baseRecv(baseCount);
		std::vector<std::array<block, 2>> baseSend(baseCount);
		BitVector baseChoice(baseCount);
		baseChoice.randomize(prng0);
		prng0.get((u8*)baseSend.data()->data(), sizeof(block) * 2 * baseSend.size());
		for (u64 i = 0; i < baseCount; ++i)
		{
			baseRecv[i] = baseSend[i][baseChoice[i]];
		}

		// set up networking
		std::string name = "n";
		IOService ios;
		Endpoint ep0(ios, "localhost", 1212, EpMode::Server, name);
		Endpoint ep1(ios, "localhost", 1212, EpMode::Client, name);
		auto recvChl = ep1.addChannel(name, name);
		auto sendChl = ep0.addChannel(name, name);


		// set the base OTs
		sender.setBaseOts(baseRecv, baseChoice);
		recv.setBaseOts(baseSend);

		u64 stepSize = 10;
		std::vector<block> inputs(stepSize);

		for (size_t j = 0; j < 2; j++)
		{
			// perform the init on each of the classes. should be performed concurrently
			auto thrd = std::thread([&]() { sender.init(numOTs, prng0, sendChl); });
			recv.init(numOTs, prng1, recvChl);
			thrd.join();

			std::vector<block> encoding1(stepSize), encoding2(stepSize);

			// Get the random OT messages
			for (u64 i = 0; i < numOTs; i += stepSize)
			{

				prng0.get(inputs.data(), inputs.size());


				for (u64 k = 0; k < stepSize; ++k)
				{

					// The receiver MUST encode before the sender. Here we are only calling encode(...) 
					// for a single i. But the receiver can also encode many i, but should only make one 
					// call to encode for any given value of i.
					recv.encode(i + k, &inputs[k], (u8*)&encoding1[k], sizeof(block));
				}

				// This call will send to the other party the next "stepSize" corrections to the sender.
				// If we had made more or less calls to encode above (for contigious i), then we should replace
				// stepSize with however many calls we made. In an extreme case, the reciever can perform
				// encode for i \in {0, ..., numOTs - 1}  and then call sendCorrection(recvChl, numOTs).
				recv.sendCorrection(recvChl, stepSize);

				// receive the next stepSize correction values. This allows the sender to now call encode
				// on the next stepSize OTs.
				sender.recvCorrection(sendChl, stepSize);

				for (u64 k = 0; k < stepSize; ++k)
				{
					// the sender can now call encode(i, ...) for k \in {0, ..., i}. 
					// Lets encode the same input and then we should expect to
					// get the same encoding.
					sender.encode(i + k, &inputs[k], (u8*)&encoding2[k], sizeof(block));

					// check that we do in fact get the same value
					if (neq(encoding1[k], encoding2[k]))
						throw UnitTestFail();

					// In addition to the sender being able to obtain the same value as the receiver,
					// the sender can encode and other codeword. This should result in a different 
					// encoding.
					inputs[k] = prng0.get<block>();

					sender.encode(i + k, &inputs[k], (u8*)&encoding2[k], sizeof(block));

					if (eq(encoding1[k], encoding2[k]))
						throw UnitTestFail();
				}
			}
		}

		// Double check that we can call split and perform the same tests.
		auto recv2Ptr = recv.split();
		auto send2Ptr = sender.split();

		auto& recv2 = *recv2Ptr;
		auto& send2 = *send2Ptr;

		for (size_t j = 0; j < 2; j++)
		{
			auto thrd = std::thread([&]() {
				send2.init(numOTs, prng0, sendChl);
			});

			recv2.init(numOTs, prng1, recvChl);

			thrd.join();


			for (u64 i = 0; i < numOTs; ++i)
			{
				block input = prng0.get<block>();

				block encoding1, encoding2;
				recv2.encode(i, &input, &encoding1);

				recv2.sendCorrection(recvChl, 1);
				send2.recvCorrection(sendChl, 1);

				send2.encode(i, &input, &encoding2);

				if (neq(encoding1, encoding2))
					throw UnitTestFail();

				input = prng0.get<block>();

				send2.encode(i, &input, &encoding2);

				if (eq(encoding1, encoding2))
					throw UnitTestFail();
			}

		}

		sendChl.close();
		recvChl.close();

		ep0.stop();
		ep1.stop();
		ios.stop();
	}




}