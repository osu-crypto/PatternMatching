#include "stdafx.h"
#ifdef  _MSC_VER
#include "CppUnitTest.h"
#include "PM_Tests.h"
#include "NcoOT_Tests.h"
#include "Common.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace tests_libOTe
{
    TEST_CLASS(OT_Tests)
    {
    public:

        TEST_METHOD(Transpose_TestVS)
        {
            InitDebugPrinting();
            PM_Transpose_Test_Impl();
        }

        TEST_METHOD(TransposeMatrixView_TestVS)
        {
            InitDebugPrinting();
			PM_TransposeMatrixView_Test_Impl();
        }

        TEST_METHOD(Iknp_200Receive_TestVS)
        {
            InitDebugPrinting();
			PM_IknpOtExt_100Receive_Test_Impl();
        }


        TEST_METHOD(Kkrt_200Receive_TestVS)
        {
            InitDebugPrinting();
            KkrtNcoOt_Test_Impl();
        }
		TEST_METHOD(MP_Matrix_Test)
		{
			InitDebugPrinting();
			MP_Matrix_Test_Impl();
		}

		
		TEST_METHOD(Shift_Test)
		{
			InitDebugPrinting();
			Shift_Test_Impl();
		}
		
		

    };
}
#endif