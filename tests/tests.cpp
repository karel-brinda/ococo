#include "ococo.h"
#include "gtest/gtest.h"


vector<int> nucls{nt16_A,nt16_C,nt16_G,nt16_T};

namespace {
	
	// The fixture for testing class Foo.
	class CounterTest : public ::testing::Test {
	protected:
		// You can remove any or all of the following functions if its body
		// is empty.
		
		CounterTest() {
			// You can do set-up work for each test here.
		}
		
		virtual ~CounterTest() {
			// You can do clean-up work that doesn't throw exceptions here.
		}
		
		// If the constructor and destructor are not enough for setting up
		// and cleaning up each test, you can define the following methods:
		
		virtual void SetUp() {
			// Code here will be called immediately after the constructor (right
			// before each test).
		}
		
		virtual void TearDown() {
			// Code here will be called immediately after each test (right
			// before the destructor).
		}
		
		// Objects declared here can be used by all tests in the test case for Foo.
	};
	
	TEST_F(CounterTest, AllIncrements) {
		counter_t c1;

		for(int nucl : nucls)
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
			c1 = _COUNTER_CELL_SET(c1,nt16_A,cell_maxval);
			c1 = _COUNTER_CELL_SET(c1,nt16_C,cell_maxval);
			c1 = _COUNTER_CELL_SET(c1,nt16_G,cell_maxval);
			c1 = _COUNTER_CELL_SET(c1,nt16_T,cell_maxval);

			c1 = _COUNTER_CELL_INC(c1,nucl);

			// c2 := max values after overflowing, then pivot
			c2=0;
			c2 = _COUNTER_CELL_SET(c1,nt16_A,cell_maxval_shifted);
			c2 = _COUNTER_CELL_SET(c1,nt16_C,cell_maxval_shifted);
			c2 = _COUNTER_CELL_SET(c1,nt16_G,cell_maxval_shifted);
			c2 = _COUNTER_CELL_SET(c1,nt16_T,cell_maxval_shifted);

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
	
	
}  // namespace
	


int main(int argc, char** argv){
	cout << "DEBUG INFO" << endl;

	cout << "cell_maxval 0x" << hex << cell_maxval << endl;
	cout << "cell_maxval_shifted 0x" << hex  << cell_maxval_shifted << endl;
	cout << "counter_norm_mask 0x" << hex  << counter_norm_mask << endl << endl;
 
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}