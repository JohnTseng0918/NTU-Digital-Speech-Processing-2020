#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <limits>
#include <algorithm>
#include "Ngram.h"

using namespace std;

double getBigramProb(const char *w1, const char *w2, Vocab& voc, Ngram& lm);

vector<string> viterbi(vector<string> textvec, map<string, vector<string>>& map_table, Vocab& voc, Ngram& lm){
	vector<vector<double>> prob_table;
	vector<vector<string>> str_table;
	vector<vector<int>> route_table;

	// init
	vector<double> column;
	vector<string> s = map_table[textvec[0]];
	vector<int> route;
	column.resize(s.size());
	route.resize(s.size());
	for (int i = 0; i < s.size(); i++){
		column[i] = getBigramProb("", s[i].c_str(), voc, lm);
		route[i] = 0;
	}
	prob_table.push_back(column);
	str_table.push_back(s);
	route_table.push_back(route);

	// induction
	
	for (int i = 1; i < textvec.size(); i++){
		column.clear();
		s = map_table[textvec[i]];
		route.clear();
		column.resize(s.size());
		route.resize(s.size());

		for (int j = 0; j < s.size(); j++){
			double max_prob = numeric_limits<double>::lowest();
			double p;
			int r;
			for (int k = 0; k < str_table[i-1].size(); k++){
				p = getBigramProb(str_table[i-1][k].c_str(), s[j].c_str(), voc, lm) + prob_table[i-1][k];
				if (p > max_prob){
					max_prob = p;
					r = k;
				}
			}
			column[j] = max_prob;
			route[j] = r;
		}

		prob_table.push_back(column);
		str_table.push_back(s);
		route_table.push_back(route);
	}

	// the last word and find index
	int idx;
	double max_prob = numeric_limits<double>::lowest();
	for (int i = 0; i < str_table.back().size(); i++){
		double p = getBigramProb(str_table.back()[i].c_str(), "", voc, lm) + prob_table.back()[i];
		if (p > max_prob){
			p = max_prob;
			idx = i;
		}
	}
	
	// trace
	vector<string> textdone;
	for (int i = str_table.size() - 1; i >= 0; i--){
		textdone.push_back(str_table[i][idx]);
		idx = route_table[i][idx];		
	}

	reverse(textdone.begin(), textdone.end());
	
	return textdone;
}


int main(int argc, char *argv[]){

	// check argv
	if (argc != 5){
		cout << "there must be 5 arguments" << endl;
		return 0;
	}

	
	char *seg_filename = argv[1];
	char *zbmap_filename = argv[2];
	char *lm_filename = argv[3];
	char *output_filename = argv[4];

	// use library load language model
	Vocab voc;
	Ngram lm(voc, 2);
	File lmFile( lm_filename, "r" );
	lm.read(lmFile);
	lmFile.close();
    
	// read ZhuYin-Big5 mapping
	map<string, vector<string>> map_table;
	ifstream mapfile(zbmap_filename, ifstream::in);

	string line;
	while(getline(mapfile, line)){
		string word = line.substr(0, 2);
		vector<string> list;
		for (int i = 3; i < line.length(); i+=2){
			list.push_back(line.substr(i, 2));
			i++;
		}
		map_table[word] = list;
	}
	mapfile.close();
	
	// read segment file
	vector<vector<string>> seg_text;
	ifstream segfile(seg_filename,ifstream::in);

	while(getline(segfile, line)){
		string tmptext = "";
		for (int i = 0; i < line.length(); i++){
			if(line[i] != ' '){
				tmptext += line[i];
			}
		}

		vector<string> tmpV;
		for (int i = 0; i < tmptext.length(); i += 2){
			tmpV.push_back(tmptext.substr(i, 2));
		}

		seg_text.push_back(tmpV);
	}
	segfile.close();

	// process segment file and output

	ofstream outputfile(output_filename, ofstream::out);

	for (vector<vector<string>>::iterator it = seg_text.begin(); it != seg_text.end(); it++){
		vector<string> textvec = *it;
		
		/*
		outputfile << "<s> ";
		for (vector<string>::iterator iter = textvec.begin(); iter!=textvec.end();iter++){
			outputfile << *iter << " ";
		}
		outputfile << "</s>" << endl;
		*/

		vector<string> textdone = viterbi(textvec, map_table, voc, lm);
		outputfile << "<s> ";
		for (vector<string>::iterator iter = textdone.begin(); iter != textdone.end(); iter++){
			outputfile << *iter << " ";
		}
		outputfile << "</s>" << endl;
	}

	outputfile.close();

	return 0;
}

// Get P(W2 | W1) -- bigram
double getBigramProb(const char *w1, const char *w2, Vocab& voc, Ngram& lm){
	VocabIndex wid1 = voc.getIndex(w1);
	VocabIndex wid2 = voc.getIndex(w2);

	//cout << wid1 << " " << wid2 << endl;

	if(wid1 == Vocab_None)  //OOV
		wid1 = voc.getIndex(Vocab_Unknown);
	if(wid2 == Vocab_None)  //OOV
		wid2 = voc.getIndex(Vocab_Unknown);
	
	VocabIndex context[] = { wid1, Vocab_None };
	return lm.wordProb( wid2, context);
}
