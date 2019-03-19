#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <cmath>
#include <utility>
#include <iterator>
#include <fstream>
#include <ctime>

const double g = 9.81;

void replacePath(std::vector <double> &coord, double coordFinal) {
	int size = coord.size();
	for (int i = 0; i < size; i++) {
		coord[i] -= coordFinal;
	}
	return;
}

std::vector <std::pair <double, double>> searchTimeSector(std::map <double, double> W, std::vector <double> y, const double h) {
	std::map <double, double> ::reverse_iterator itrStart = W.rbegin();
	std::map <double, double> ::reverse_iterator itrEnd = --W.rend();
	std::map <double, double> ::iterator itrY1, itrY2;
	std::vector <std::pair <double, double>> time;
	double stpH, currH = h;
	for (std::map <double, double> ::reverse_iterator it = itrStart; it != itrEnd; ++it) {
		stpH = it->first - (++it)->first;
		--it;
		int j1 = 0;
		while (y[j1] > currH) {
			j1++;
		}
		int j2 = j1;
		if (currH - stpH > 0) {
			while (y[j2] > (currH - stpH)) {
				j2++;
			}
		}
		else {
			j2 = y.size();
		}
		time.push_back(std::pair <double, double>(currH - stpH, j2 - j1));
		currH -= stpH;
	}
	return time;
}

std::pair <double, double> interpolateForce(std::map <double, double> aeroPower, const double speed) {
	std::map <double, double> ::iterator itr = aeroPower.lower_bound(speed);
	double y1 = (itr)->second;
	double x1 = (itr)->first;
	double y2 = (++itr)->second;
	double x2 = (itr)->first;
	double k = (y2 - y1) / (x2 - x1);
	double b = (x2*y1 - x1*y2) / (x2 - x1);
	return std::pair <double, double>(k, b);
}

std::vector<double> yCoord(const double m, const double dt, std::map <double, double> aeroPower, double h, std::vector<double> &ySpeed) {
	std::vector <double> y;
	double k, b, a, speed = 0;
	std::map <double, double> ::iterator itr = ++aeroPower.begin();
	std::pair <double, double> lineForce;
	while (h > 0)
	{
		y.push_back(h);
		ySpeed.push_back(speed);
		if (speed > itr->first) {
			lineForce = interpolateForce(aeroPower, (speed));
			itr++;
		}
		a = abs(-g*m + lineForce.first*speed + lineForce.second) / m;
		h -= speed*dt + a*dt*dt / 2;
		speed = speed + a * dt;
	}
	return (y);
}

std::vector <double> coordinates(const double m, const double dt, std::map <double, double> Wx, std::map <double, double> aeroPower, double xSpeed,
	std::vector <double> y, std::vector <std::pair<double, double>> timeSector, std::vector <double> &coordSpeed) {
	std::vector <double> x;
	double k, b, aBody, currX = 0, aWind;
	int number1 = timeSector[0].second, number2 = 0, number3 = 0; //!!!!!!!!!!!!!!!
	std::map <double, double> ::iterator itr = ++aeroPower.begin();
	int y_size = y.size();
	std::pair <double, double> lineForce;
	for (int i = 0; i < y_size; i++)
	{
		x.push_back(currX);
		coordSpeed.push_back(xSpeed);
		if (xSpeed > itr->first) {
			lineForce = interpolateForce(aeroPower, (xSpeed));
			itr++;
		}
		std::map <double, double> ::iterator itr_Wx = Wx.lower_bound(y[i]);
		std::vector <std::pair<double, double>> ::iterator itrTimeSector = timeSector.begin();
		/*while (itrTimeSector->first > y[i]) {
			if (y[i] > 200) {
				++itrTimeSector;
			}
			else {
				itrTimeSector = --timeSector.end();
				break;
			}
		}
		aWind = (-itr_Wx->second + (--itr_Wx)->second) / (dt*itrTimeSector->second);*/
		int k = 0;
		for (int j = 0; j<y_size, timeSector[j].first > y[i]; j++) {
			k++;
		}
		aWind = (-itr_Wx->second + (--itr_Wx)->second) / (dt*timeSector[k].second);

		/*if (number1 != itrTimeSector ->second) {
			number1 = itrTimeSector->second;
			std::cout << number2 << std::endl;
			number2 = 0;
		}
		number2++;*/

		aBody = (aWind*m - (lineForce.first*xSpeed + lineForce.second)) / m;
		currX += xSpeed*dt + aBody*dt*dt / 2;
		xSpeed = xSpeed + aBody * dt;
		//std::cout <<"time " << timeSector[k].second << ' ' <<currX<<' '<< x_speed<<" aWind:" <<aWind <<" forse:"<<lineForce.first*x_speed + lineForce.second <<" aBody:"<<aBody<<std::endl;
	}
	return (x);
}

void read_F_csv(std::map<double, double> &F_map, std::string file_name) {
	std::ifstream f(file_name);
	char dummy;
	float a, b;
	while (f >> a >> dummy >> b)
	{
		F_map[a] = b;
	}
	f.close();
	return;
}

void read_Wind_csv(std::map<double, double> &X_map, std::map<double, double> &Z_map, std::string file_name) {
	std::ifstream f(file_name);
	char dummy1, dummy2;
	float Y, Wx, Wz;
	while (f >> Y >> dummy1 >> Wx >> dummy2 >> Wz)
	{
		X_map[Y] = Wx;
		Z_map[Y] = Wz;
	}
	f.close();
	return;
}

void CSV_output(std::map<float, float> F_map) {
	for (auto it = F_map.begin(); it != F_map.end(); ++it)
	{
		std::cout << (*it).first << " : " << (*it).second << std::endl;
	}
	system("pause");
	system("cls");
}

int main() {
	unsigned int start_time = clock();
	std::cout.precision(10);
	double dt = 0.01, h = 1400, m = 100, startSpeed = 250, dfi = 2 * M_PI / 360, speed;
	//std::cin >> "¬ведите скорость" >> speed >> "¬ведите высоту" >> h;

	std::map <double, double> aeroPower, Wx, Wz;
	read_F_csv(aeroPower, "F.csv");
	read_Wind_csv(Wx, Wz, "Wind.csv");
	std::vector <double> x, z;
	std::vector <double>  y_speed, x_speed, z_speed;;
	std::vector <double> y = yCoord(m, dt, aeroPower, h, y_speed);
	int size_of_coord = y.size();
	std::vector <std::pair<double, double>> time = searchTimeSector(Wx, y, h);
	double  M[89][1]; 
	double alf=0;
	for (int i = 0; i < 90; i++) {
		if (i == 89) {
		//	system("pause");
		}
		alf = alf + 2 * M_PI / 360;
		x = coordinates(m, dt, Wx, aeroPower, startSpeed*cos(alf), y, time, x_speed);
		z = coordinates(m, dt, Wz, aeroPower, startSpeed*sin(alf), y, time, z_speed);
		x_speed.clear();
		z_speed.clear();
		std::cout <<alf<<' '  <<startSpeed*cos(alf)<<' ' << x[x.size() - 1] << ' ' << z[z.size() - 1] << ' ' << std::endl;
		replacePath(x, x[x.size() - 1]);
		replacePath(z, z[z.size() - 1]);
		M[i][0] = x[0]; M[i][1] = z[0];
		std::cout << M[i][0] << ' ' << M[i][1] << ' ' << std::endl;
		std::cout << i << std::endl;
	}

	std::ofstream out;
	out.open("\trajectory.txt");
	if (out.is_open())
	{
		for (int i = 0; i < size_of_coord; i++) {
		//	out << dt*i << ' ' << sqrt(x_speed[i] * x_speed[i] + y_speed[i] * y_speed[i] + z_speed[i] * z_speed[i]) << ' ' << x[i] << ' ' << z[i] << ' ' << y[i] << std::endl;
		}
	}
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	std::cout << search_time << std::endl;

	return 0;
}


/*«дравствуйте ƒаниил. ¬ представленом задании на BestHack требуетс€ определить угол и начальную кординату сброса, однако, учитыва€ представленные данные (аэродинамическую силу
,скорости ветра не разных высотах, высота сброса, начальна€ скорость и массу груза), определить параметры однозначно невозможно, так как в зависимости от угла измен€етс€ и
траектори€ полета, а следовательно получаетс€ множество решений*/