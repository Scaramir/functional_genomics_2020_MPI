#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
//using namespace std;

//split a row and push the row-elements in a vector
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

//get ascii-score of the string, to check, if it's below the score of non-autosomal genes
bool autosom(std::string &a){
    bool b = false; 
    if (a.size() <= 4){
        for (int i = 0; i < a.size()-2; ++i){
            if (isdigit(a[i+1]))
                b = true;
            else 
                b = false;
        }
        return b;
    }
    return false;
}

///#####MAIN#####
int main(){
    //Input-Pfade anpassen!
    std::string linec;
    std::ifstream counts ("lung_vs_forebrain.txt");
    std::string lineg;
    std::ifstream prot_gene ("prot.gene.txt");
    //output (filtered):
    std::ofstream out;
    out.open("lung_vs_forebrain_filtered.txt");
    if (!counts.good() || !prot_gene.good()){
		std::cerr << "Unable to open requested files!" << std::endl;
		return 1; //cancel
    } else {
        std::cout << "Files can be opened \n";
        prot_gene.close();
        std::cout << "\n"; ///

        std::getline(counts, linec);
        while (std::getline(counts, linec)){
            std::vector<std::string> row_values_c;
            split(linec, '\t', row_values_c);   

            std::string lineg; //better: open file once, start getline on every while-beginning from line 1
            std::ifstream prot_gene ("prot.gene.txt");
            while (std::getline(prot_gene, lineg)){
                std::vector<std::string> row_values_g;
                split(lineg, '\t', row_values_g);
                std::string gene_id = '"' + row_values_c[0] + '"'; 
                if (gene_id == row_values_g[0]){        
                    if (autosom(row_values_g[4])){
                    out << linec << "\n";
                    }
                }
            }
        }
        prot_gene.close();
        out.close();
        std::cout << "done \n";
        return 0;
    }
    return 0;
}   	

