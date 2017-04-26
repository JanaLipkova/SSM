// vector::begin/end
#include <iostream>
#include <vector>

int main ()
{
	std::vector<int> myvector;
	for (int i=5; i<=10; i++) myvector.push_back(i*10);
	
	std::cout << "myvector contains:";
	for (std::vector<int>::iterator it = myvector.begin() ; it != myvector.end(); ++it)
		std::cout << ' ' << *it;
	std::cout << '\n';
	
	std::vector<int>::iterator it = myvector.begin();
	
	while (it != myvector.end())
	{
		std::cout <<"it"<<*it<<std::endl;
		++it;
		std::cout <<"++it"<<*it<<std::endl;
	}
	
	return 0;
}
