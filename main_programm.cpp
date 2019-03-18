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

void replace_path(std::vector <double> &coord, double coord_final) {
	int size_of_x = coord.size();
	for (int i = 0; i < size_of_x; i++) {
		coord[i] -= coord_final;
	}
	return;
}

std::vector <std::pair <double, double>> time_sector(std::map <double, double> W, std::vector <double> y, const double h) {
	std::map <double, double> ::reverse_iterator itr_start = W.rbegin();
	std::map <double, double> ::reverse_iterator itr_end = --W.rend();
	std::map <double, double> ::iterator itr_y1, itr_y2;
	std::vector <std::pair <double, double>> time;
	double stp_h, curr_h = h;
	for (std::map <double, double> ::reverse_iterator it = itr_start; it != itr_end; ++it) {
		stp_h = it->first - (++it)->first;
		--it;
		int j1 = 0;
		while (y[j1] > curr_h) {
			j1++;
		}
		int j2 = j1;
		if (curr_h - stp_h > 0) {
			while (y[j2] > (curr_h - stp_h)) {
				j2++;
			}
		}
		else {
			j2 = y.size();
		}
		time.push_back(std::pair <double, double>(curr_h - stp_h, j2 - j1));
		curr_h -= stp_h;
	}
	return time;
}

std::pair <double, double> interpolate_force(std::map <double, double> aero_power, const double speed) {
	std::map <double, double> ::iterator itr = aero_power.lower_bound(speed);
	double y1 = (itr)->second;
	double x1 = (itr)->first;
	double y2 = (++itr)->second;
	double x2 = (itr)->first;
	double k = (y2 - y1) / (x2 - x1);
	double b = (x2*y1 - x1*y2) / (x2 - x1);
	return std::pair <double, double>(k, b);
}

std::vector<double> y_coord(const double m, const double dt, std::map <double, double> aero_power, double h, std::vector<double> &y_speed) {
	std::vector <double> y;
	double k, b, a, speed = 0;
	std::map <double, double> ::iterator itr = ++aero_power.begin();
	std::pair <double, double> line_force;
	while (h > 0)
	{
		y.push_back(h);
		y_speed.push_back(speed);
		if (speed > itr->first) {
			line_force = interpolate_force(aero_power, (speed));
			itr++;
		}
		a = abs(-g*m + line_force.first*speed + line_force.second) / m;
		h -= speed*dt + a*dt*dt / 2;
		speed = speed + a * dt;
	}
	return (y);
}

std::vector <double> coordinates(const double m, const double dt, std::map <double, double> Wx, std::map <double, double> aero_power, double x_speed,
	std::vector <double> y, std::vector <std::pair<double, double>> time_sector, std::vector <double> &coord_speed) {
	std::vector <double> x;
	double k, b, a_body, curr_x = 0, a_wind;
	int number1 = time_sector[0].second, number2 = 0, number3 = 0; //!!!!!!!!!!!!!!!
	std::map <double, double> ::iterator itr = ++aero_power.begin();
	int y_size = y.size();
	std::pair <double, double> line_force;
	for (int i = 0; i < y_size; i++)
	{
		x.push_back(curr_x);
		coord_speed.push_back(x_speed);
		if (x_speed > itr->first) {
			line_force = interpolate_force(aero_power, (x_speed));
			itr++;
		}
		std::map <double, double> ::iterator itr_Wx = Wx.lower_bound(y[i]);
		std::vector <std::pair<double, double>> ::iterator itr_time_sector = time_sector.begin();
		/*while (itr_time_sector->first > y[i]) {
			if (y[i] > 200) {
				++itr_time_sector;
			}
			else {
				itr_time_sector = --time_sector.end();
				break;
			}
		}
		a_wind = (-itr_Wx->second + (--itr_Wx)->second) / (dt*itr_time_sector->second);*/
		int k = 0;
		for (int j = 0; j<y_size, time_sector[j].first > y[i]; j++) {
			k++;
		}
		a_wind = abs(-itr_Wx->second + (--itr_Wx)->second) / (dt*time_sector[k].second);

		/*if (number1 != itr_time_sector ->second) {
			number1 = itr_time_sector->second;
			std::cout << number2 << std::endl;
			number2 = 0;
		}
		number2++;*/

		a_body = (a_wind*m - (line_force.first*x_speed + line_force.second)) / m;
		curr_x += x_speed*dt + a_body*dt*dt / 2;
		x_speed = x_speed + a_body * dt;
		//std::cout <<"time " << time_sector[k].second << ' ' <<curr_x<<' '<< x_speed<<" a_wind:" <<a_wind <<" forse:"<<line_force.first*x_speed + line_force.second <<" a_body:"<<a_body<<std::endl;
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
	double dt = 0.01, h = 1400, m = 100, start_speed = 250, dfi = 2 * M_PI / 360, speed;
	//std::cin >> "¬ведите скорость" >> speed >> "¬ведите высоту" >> h;

	std::map <double, double> aero_power, Wx, Wz;
	read_F_csv(aero_power, "F.csv");
	read_Wind_csv(Wx, Wz, "Wind.csv");
	std::vector <double> x, z, x_speed, y_speed, z_speed;
	std::vector <double> y = y_coord(m, dt, aero_power, h, y_speed);
	int size_of_coord = y.size();
	std::vector <std::pair<double, double>> time = time_sector(Wx, y, h);
	double M[90][3];
	int i = 0;
	double alf = 45;
	//for (int alf = 0; alf < M_PI / 2; alf += 2 * M_PI / 360) {
		x = coordinates(m, dt, Wx, aero_power, start_speed*cos(alf), y, time, x_speed);
		z = coordinates(m, dt, Wz, aero_power, start_speed*sin(alf), y, time, z_speed);
		/*x_speed.clear();
		z_speed.clear();*/
		replace_path(x, x[x.size() - 1]);
		replace_path(z, z[z.size() - 1]);
		M[i][0] = x[x.size() - 1]; M[i][1] = z[z.size() - 1]; M[i][2] = y[y.size() - 1];
		std::cout << i << std::endl;
		i++;
	//}

	std::ofstream out;
	out.open("\trajectory.txt");
	if (out.is_open())
	{
		for (int i = 0; i < size_of_coord; i++) {
			out << dt*i << ' ' << sqrt(x_speed[i] * x_speed[i] + y_speed[i] * y_speed[i] + z_speed[i] * z_speed[i]) << ' ' << x[i] << ' ' << z[i] << ' ' << y[i] << std::endl;
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