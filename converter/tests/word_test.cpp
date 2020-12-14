#include "test.h"
#include <vector>
#include <fstream>

int main(){
    //int test(){

    l1tk_word_t in_tk(0);
    pf_tk_word_t out_tk(0);

    std::vector<std::string> in_words{
        // add words here for testing a few examples
        // "e281083934002762fd0332bf",
    };

    // read 96b tracks from an input file
    std::ifstream infile("../../../../tests/word_test_input.txt");
    std::string word;
    while (infile >> word){
        if (word.size()!=24){ // 96b = 24 hex chars * 4b / hex char
            cout << word << " has too few characters for a 96b word. skipping!!" << endl;
            continue;
        }
        in_words.push_back( word );
    }

    // write 64b tracks to an output file
    std::ofstream outfile;
    outfile.open("../../../../tests/word_test_output.txt");

    for(auto in_word : in_words){
        l1tk_word_t in_tk(in_word.c_str(), 16);
        pf_input_track_conv_hw(in_tk, out_tk, 0/*linkNo*/);
        outfile << out_tk.to_string(16) << endl;
    }

    return 0;
}
