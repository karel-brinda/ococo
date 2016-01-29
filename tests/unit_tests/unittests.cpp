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
	
	class CounterTest : public ::testing::Test {
	protected:
		CounterTest() {
		}
		
		virtual ~CounterTest() {
		}
		
		virtual void SetUp() {
		}
		
		virtual void TearDown() {
		}
	};
	
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
	
	TEST_F(CounterTest, AllIncrements) {
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
	}

	TEST_F(NuclGeneratorTest, IndividualNucleotides) {
		char nucl;

        counters_quadruplet_t quadruplet;
        
        quadruplet={0,0,0,0,0};
		nucl=rand_nucl(quadruplet);
		ASSERT_EQ(nucl, 'N');

        quadruplet={5,0,0,0,5};
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
		ASSERT_EQ(nucl, 'T');

	}
	
	
}  // namespace
	


int main(int argc, char** argv){
	cout << endl << "DEBUG INFO" << endl;

	cout << "\tcounter_size " << dec << counter_size << "B" << endl;
	cout << "\tcell_bits " << dec << cell_bits << endl;
	cout << "\tcell_maxval 0x" << hex << cell_maxval << endl;
	cout << "\tcell_maxval_shifted 0x" << hex  << cell_maxval_shifted << endl;
	cout << "\tcounter_norm_mask 0x" << hex  << counter_norm_mask << endl << endl;
 
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
