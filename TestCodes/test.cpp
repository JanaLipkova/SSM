 #include <iostream>     // std::cout
 #include <algorithm>    // std::find
 #include <vector>       // std::vector


using namespace std;

int main()
{
	vector<int> v;
	
	if (v.empty())
	{
		cout << "is empty" <<endl;
	}
	
	v.push_back(1);
	
	if (v.empty()) {
		cout<<" is still empty"<<endl;
	}else {
		cout <<" now is not empty"<<endl;
	}

	
	return 0;
}