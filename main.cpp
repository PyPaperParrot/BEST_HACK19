#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
//#include <algorithm>
#include <map>


//Read F.csv file into map
bool Read_F_csv(std::map<float, float> &F_map, std::string file_name) {
	
	std::ifstream f(file_name);

	char dummy;
	float a, b;
	
	while (f >> a >> dummy >> b)
	{
		F_map[a] = b;
	}

	f.close();
	return true;
}
//Read Wind.csv file into mao
bool Read_Wind_csv(std::map<float, float> &X_map, 
	std::map<float, float> &Z_map, std::string file_name) {

	std::ifstream f(file_name);

	char dummy1, dummy2;
	float Y, Wx, Wz;

	while(f >> Y >> dummy1 >> Wx >> dummy2 >> Wz)
	{
		X_map[Y] = Wx;
		Z_map[Y] = Wz;
	}

	f.close();
	return true;
}
//output maps
void CSV_output(std::map<float, float> F_map) {
	for (auto it = F_map.begin(); it != F_map.end(); ++it)
	{
		std::cout << (*it).first << " : " << (*it).second << std::endl;
	}
	system("pause");
	system("cls");
}


int main() {
	
	std::map<float, float> F_map;
	std::map<float, float> Wx_map;
	std::map<float, float> Wz_map;
	
	Read_F_csv(F_map, "F.csv");
	Read_Wind_csv(Wx_map, Wz_map, "Wind.csv");

	CSV_output(F_map);
	CSV_output(Wx_map);
	CSV_output(Wz_map);
	
	return 0;
}