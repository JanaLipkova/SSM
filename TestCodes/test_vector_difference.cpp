// set_difference example
#include <iostream>     // std::cout
#include <algorithm>    // std::set_difference, std::sort
#include <vector>       // std::vector

using namespace std;

int main () {
	
	vector<int> all;
	vector<int> crit;
	
	for(int i=0; i<5; i++)
		all.push_back(i);
	
	crit.push_back(1);
	crit.push_back(3);
	

	std::vector<int> v(10);                      // 0  0  0  0  0  0  0  0  0  0
	std::vector<int>::iterator it;

	
	it=std::set_difference (all.begin(), all.end(), crit.begin(), crit.end(), v.begin());
// 5 15 25  0  0  0  0  0  0  0
	v.resize(it-v.begin());                      //  5 15 25

	std::cout << "The difference has " << (v.size()) << " elements:\n";
	for (it=v.begin(); it!=v.end(); ++it)
		std::cout << ' ' << *it;
	std::cout << '\n';
	
	return 0;
}