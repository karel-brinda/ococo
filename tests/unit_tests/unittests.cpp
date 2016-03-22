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
    
	
	class ConsensusTest : public ::testing::Test {
	protected:
		ConsensusTest() {
		}
		
		virtual ~ConsensusTest() {
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
       
    }
    
	TEST_F(ConsensusTest, EmptyStats) {
		char nucl;

        params_t params=params_t();
        for(uint32_t i=0;i<4;i++){
            {
                ococo::pos_stats_uncompr_t psu = {nt256_nt16[(int)'C'],{0,0,0,0},0};
                nucl=(params.cons_alg[i])(psu, params);
                ASSERT_EQ('C',nucl);
            }
            
            {
                ococo::pos_stats_uncompr_t psu = {nt256_nt16[(int)'N'],{0,0,0,0},0};
                nucl=(params.cons_alg[i])(psu, params);
                ASSERT_EQ('N',nucl);
            }
        }
        
	}
	
    TEST_F(ConsensusTest, Majority) {
        char nucl;
        
        params_t params=params_t();
        params.majority_threshold=0.6;

        {
            ococo::pos_stats_uncompr_t psu = {nt256_nt16[(int)'T'],{6,0,0,3},9};
            nucl=cons_call_maj(psu, params);
            ASSERT_EQ('A',nucl);
        }
        
        {
            ococo::pos_stats_uncompr_t psu = {nt256_nt16[(int)'T'],{0,0,6,4},10};
            nucl=cons_call_maj(psu, params);
            ASSERT_EQ('G',nucl);
        }
        
    }
    
	
}  // namespace
	


int main(int argc, char** argv){
	cout << endl << "DEBUG INFO" << endl;

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
