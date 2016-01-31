#include "ococo.h"
#include "gtest/gtest.h"

#include <vector>

using namespace ococo;
using namespace std;

vector<uint32_t> nucls {
    nt256_nt16[(int)'A'],
    nt256_nt16[(int)'C'],
    nt256_nt16[(int)'G'],
    nt256_nt16[(int)'T']
};

namespace {

    class BitFunctionsTest : public ::testing::Test {
    protected:
        BitFunctionsTest() {
        }
        
        virtual ~BitFunctionsTest() {
        }
        
        virtual void SetUp() {
        }
        
        virtual void TearDown() {
        }
    };
    
	/*class CounterTest : public ::testing::Test {
	protected:
		CounterTest() {
		}
		
		virtual ~CounterTest() {
		}
		
		virtual void SetUp() {
		}
		
		virtual void TearDown() {
		}
	};*/
	
	class NuclGeneratorTest : public ::testing::Test {
	protected:
		NuclGeneratorTest() {
		}
		
		virtual ~NuclGeneratorTest() {
		}
		
		virtual void SetUp() {
		}
		
		virtual void TearDown() {
		}
	};

    TEST_F(BitFunctionsTest, Basic) {
        ASSERT_EQ( 0x00,   (ococo::right_full_mask<uint8_t,0>()) );
        ASSERT_EQ( 0x01,   (ococo::right_full_mask<uint8_t,1>()) );
        ASSERT_EQ( 0xff,   (ococo::right_full_mask<uint8_t,8>()) );
        ASSERT_EQ( 0xffff, (ococo::right_full_mask<uint16_t,16>()) );
        
        ASSERT_EQ( 0x00,   (ococo::left_full_mask<uint8_t,0>()) );
        ASSERT_EQ( 0x80,   (ococo::left_full_mask<uint8_t,1>()) );
        ASSERT_EQ( 0xff,   (ococo::left_full_mask<uint8_t,8>()) );
        ASSERT_EQ( 0xffff, (ococo::left_full_mask<uint16_t,16>()) );
        
        
        ASSERT_EQ( 0x01,   (ococo::get_right_bits<uint8_t,1,0>(0xff)) );
        ASSERT_EQ( 0x01,   (ococo::get_right_bits<uint8_t,1,1>(0xff)) );
        ASSERT_EQ( 0x03,   (ococo::get_right_bits<uint8_t,2,0>(0xff)) );
        ASSERT_EQ( 0xff,   (ococo::get_right_bits<uint8_t,8,0>(0xff)) );
        ASSERT_EQ( 0x08,   (ococo::get_right_bits<uint8_t,4,0>(0xf8)) );

        ASSERT_EQ( 0x01,   (ococo::get_left_bits<uint8_t,1,0>(0xff)) );
        ASSERT_EQ( 0x01,   (ococo::get_left_bits<uint8_t,1,1>(0xff)) );
        ASSERT_EQ( 0x03,   (ococo::get_left_bits<uint8_t,2,0>(0xff)) );
        ASSERT_EQ( 0xff,   (ococo::get_left_bits<uint8_t,8,0>(0xff)) );
    }
    
    
    
	/*TEST_F(CounterTest, AllIncrements) {
		counter_t c1;

		for(uint32_t nucl : nucls)
		{
			c1 =0;
			for(int i=1;i<=cell_maxval;i++){
				c1 = _COUNTER_CELL_INC(c1,nucl);
				ASSERT_EQ( _COUNTER_CELL_VAL(c1,nucl),i );
			}
		}
	}
	
	TEST_F(CounterTest, Normalization_MaxVsMax) {
		counter_t c1;
		counter_t c2;

		for(int nucl : nucls)
		{
			// c1 := max possible values, then overflow
			c1=0;
			c1 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'A'],cell_maxval);
			c1 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'C'],cell_maxval);
			c1 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'G'],cell_maxval);
			c1 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'T'],cell_maxval);

			c1 = _COUNTER_CELL_INC(c1,nucl);

			// c2 := max values after overflowing, then pivot
			c2=0;
			c2 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'A'],cell_maxval_shifted);
			c2 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'C'],cell_maxval_shifted);
			c2 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'G'],cell_maxval_shifted);
			c2 = _COUNTER_CELL_SET(c1,nt256_nt16[(int)'T'],cell_maxval_shifted);

			c2 = _COUNTER_CELL_SET(c2,nucl,1<<(cell_bits-1));


			ASSERT_EQ(c1, c2);
		}
	}

	TEST_F(CounterTest, Normalization_MaxVsOnes) {
		counter_t c1;
		counter_t c2;

		for(int nucl : nucls)
		{
			c1=0;
			// 1 to everything
			for(int nucl1 : nucls)
			{
				c1 = _COUNTER_CELL_INC(c1,nucl1);
			}
			// max value
			c1 = _COUNTER_CELL_SET(c1,nucl,cell_maxval);
			// overflow
			c1 = _COUNTER_CELL_INC(c1,nucl);
			c2 = _COUNTER_CELL_SET(0,nucl,1<<(cell_bits-1));

			ASSERT_EQ(c1, c2);
		}
	}*/

	TEST_F(NuclGeneratorTest, IndividualNucleotides) {
		char nucl;

        ococo::pos_stats_uncompr_t psu = {'A',{0,0,0,0},0};
        nucl=cons_call_stoch(psu);
        ASSERT_EQ(nucl, 'N');
        
        /*quadruplet={5,0,0,0,5};
        nucl=rand_nucl(quadruplet);
		ASSERT_EQ(nucl, 'A');

        quadruplet={0,5,0,0,5};
        nucl=rand_nucl(quadruplet);
		ASSERT_EQ(nucl, 'C');

        quadruplet={0,0,5,0,5};
        nucl=rand_nucl(quadruplet);
		ASSERT_EQ(nucl, 'G');

        quadruplet={0,0,0,5,5};
        nucl=rand_nucl(quadruplet);
		ASSERT_EQ(nucl, 'T');*/

	}
	
	
}  // namespace
	


int main(int argc, char** argv){
	cout << endl << "DEBUG INFO" << endl;

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
