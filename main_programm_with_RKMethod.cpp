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

void replacePath(std::vector <double> &coordinates, double coordinatesFinal) {
	int size = coordinates.size();
	for (int i = 0; i < size; i++) {
		coordinates[i] -= coordinatesFinal;
	}
	return;
}

std::vector <std::pair <double, double>> searchTimeSector(std::map <double, double> W, std::vector <double> y, const double h) {
	std::map <double, double> ::reverse_iterator itrStart = W.rbegin();
	std::map <double, double> ::reverse_iterator itrEnd = --W.rend();
	std::map <double, double> ::iterator itrY1, itrY2;
	std::vector <std::pair <double, double>> time;
	double stpH, currentH = h;
	for (std::map <double, double> ::reverse_iterator it = itrStart; it != itrEnd; ++it) {
		stpH = it->first - (++it)->first;
		--it;
		int j1 = 0;
		while (y[j1] > currentH) {
			j1++;
		}
		int j2 = j1;
		if (currentH - stpH > 0) {
			while (y[j2] > (currentH - stpH)) {
				j2++;
			}
		}
		else {
			j2 = y.size();
		}
		time.push_back(std::pair <double, double>(currentH - stpH, j2 - j1));
		currentH -= stpH;
	}
	return time;
}

std::pair <double, double> interpolateForce(std::map <double, double> aeroForce, const double speed) {
	std::map <double, double> ::iterator itr = aeroForce.lower_bound(speed);
	double y1 = (itr)->second;
	double x1 = (itr)->first;
	double y2 = (++itr)->second;
	double x2 = (itr)->first;
	double k = (y2 - y1) / (x2 - x1);
	double b = (x2*y1 - x1*y2) / (x2 - x1);
	return std::pair <double, double>(k, b);
}

std::vector<double> computeYCoordEulerMethod(const double m, const double dt, std::map <double, double> aeroForce, double h, std::vector<double> &ySpeed) {
	std::vector <double> y;
	double k, b, a, speed = 0;
	std::map <double, double> ::iterator itr = ++aeroForce.begin();
	std::pair <double, double> lineForce;
	while (h > 0)
	{
		y.push_back(h);
		ySpeed.push_back(speed);
		if (abs(speed) > itr->first) {
			lineForce = interpolateForce(aeroForce, (abs(speed)));
			itr++;
		}
		/*a = abs(-g*m + lineForce.first*speed + lineForce.second) / m;

		h -= speed*dt + a*dt*dt / 2;
		speed = speed + a * dt;*/
		a = (-g*m + lineForce.first*abs(speed) + lineForce.second) / m;
		h += speed*dt + a*dt*dt / 2;
		speed = speed + a * dt;
	}
	return (y);
}

std::vector<double> computeYCoordRKMethod(const double m, const double dt, std::map <double, double> aeroForce, double h, std::vector<double> &ySpeed) {
	std::vector <double> y;
	double k, b, a;
	double speed = 0, speedCorrection = 0;
	std::map <double, double> ::iterator itr = ++aeroForce.begin();
	std::pair <double, double> lineForce;
	while (h > 0)
	{
		y.push_back(h);
		ySpeed.push_back(speed);
		if (abs(speed) > itr->first) {
			lineForce = interpolateForce(aeroForce, abs(speed));
			itr++;
		}
		/*a = abs(-g*m + lineForce.first*speed + lineForce.second) / m;
		h -= speed*dt + a*dt*dt / 2;
		speedCorrection = speed + a * dt;
		speed = speed + dt*((abs(-g*m + lineForce.first*speed + lineForce.second)/m+ abs(-g*m + lineForce.first*speedCorrection + lineForce.second)/m)/2);*/

		a = (-g*m + lineForce.first*abs(speed) + lineForce.second) / m;
		h += speed*dt + a*dt*dt / 2;
		speedCorrection = speed + a * dt;
		speed = speed + dt*(((-g*m + lineForce.first*abs(speed) + lineForce.second) / m + (-g*m + lineForce.first*abs(speedCorrection) + lineForce.second) / m) / 2);
	}
	return (y);
}

std::vector <double> computeCoordinatesEulerMethod(const double m, const double dt, std::map <double, double> Wx, std::map <double, double> aeroForce, double xSpeed,
	std::vector <double> y, std::vector <std::pair<double, double>> timeSector, std::vector <double> &coordSpeed) {
	std::vector <double> x;
	double k, b, aBody, currentX = 0, aWind;
	std::map <double, double> ::iterator itr = ++aeroForce.begin();
	int y_size = y.size();
	std::pair <double, double> lineForce;
	for (int i = 0; i < y_size; i++)
	{
		x.push_back(currentX);
		coordSpeed.push_back(xSpeed);
		if (xSpeed > itr->first) {
			lineForce = interpolateForce(aeroForce, (xSpeed));
			itr++;
		}
		std::map <double, double> ::iterator itr_Wx = Wx.lower_bound(y[i]);
		std::vector <std::pair<double, double>> ::iterator itrTimeSector = timeSector.begin();

		int k = 0;
		for (int j = 0; j<y_size, timeSector[j].first > y[i]; j++) {
			k++;
		}
		aWind = (-itr_Wx->second + (--itr_Wx)->second) / (dt*timeSector[k].second);
		aBody = (aWind*m - (lineForce.first*xSpeed + lineForce.second)) / m;
		currentX += xSpeed*dt + aBody*dt*dt / 2;
		xSpeed = xSpeed + aBody * dt;
		//std::cout <<"time " << timeSector[k].second << ' ' <<currentX<<' '<< x_speed<<" aWind:" <<aWind <<" forse:"<<lineForce.first*x_speed + lineForce.second <<" aBody:"<<aBody<<std::endl;
	}
	return (x);
}

std::vector <double> computeCoordinatesRKMethod(const double m, const double dt, std::map <double, double> Wx, std::map <double, double> aeroForce, double xSpeed,
	std::vector <double> y, std::vector <std::pair<double, double>> timeSector, std::vector <double> &coordSpeed) {
	std::vector <double> x;
	double k, b, aBody, currentX = 0, aWind, xSpeedCorrection = xSpeed;
	std::map <double, double> ::iterator itr = ++aeroForce.begin();
	int y_size = y.size();
	std::pair <double, double> lineForce;
	for (int i = 0; i < y_size; i++)
	{
		x.push_back(currentX);
		coordSpeed.push_back(xSpeed);
		if (xSpeed > itr->first) {
			lineForce = interpolateForce(aeroForce, (xSpeed));
			itr++;
		}
		std::map <double, double> ::iterator itr_Wx = Wx.lower_bound(y[i]);

		int k = 0;
		for (int j = 0; j<y_size, timeSector[j].first > y[i]; j++) {
			k++;
		}
		aWind = (-itr_Wx->second + (--itr_Wx)->second) / (dt*timeSector[k].second);
		aBody = (aWind*m - (lineForce.first*xSpeed + lineForce.second)) / m;
		currentX += xSpeed*dt + aBody*dt*dt / 2;
		xSpeedCorrection = xSpeed + aBody * dt;
		xSpeed = xSpeed + dt*(aBody + (aWind*m - (lineForce.first*xSpeedCorrection + lineForce.second)) / m) / 2;
		//std::cout <<"time " << timeSector[k].second << ' ' <<currentX<<' '<< x_speed<<" aWind:" <<aWind <<" forse:"<<lineForce.first*x_speed + lineForce.second <<" aBody:"<<aBody<<std::endl;
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

void writeFile(const std::vector<double> vec, std::string file_name) {
	std::ofstream out;
	out.open(file_name);
	if (out.is_open())
	{
		//for (int i = 0; i < vec.size(); i++) {
				//out << dt*i << ' ' << sqrt(x_speed[i] * x_speed[i] + y_speed[i] * y_speed[i] + z_speed[i] * z_speed[i]) << ' ' << x[i] << ' ' << z[i] << ' ' << y[i] << std::endl;
			std::copy(vec.begin(), vec.end(), std::ostream_iterator<double>(out, "\n"));
		//}
	}
	out.close();
	return;
}

int main() {
	setlocale(0, "russian");
	std::cout.precision(10);

	unsigned int start_time = clock();
	double dt = 0.01;
	double h = 1400.0;
	double m = 100.0;
	double startSpeed = 250.0;
	double dfi = 2 * M_PI / 360, speed;

	bool isSingleTrajectory;
	std::cout << "Задать одиночную траекторию ? ";
	std::cin >> isSingleTrajectory;

	std::map <double, double> aeroForce, Wx, Wz;
	read_F_csv(aeroForce, "F.csv");
	read_Wind_csv(Wx, Wz, "Wind.csv");
	std::vector <double> x, z;
	std::vector <double>  y_speed, x_speed, z_speed;

	unsigned int start_time_computeY = clock();
	std::vector <double> y = computeYCoordRKMethod(m, dt, aeroForce, h, y_speed);
	std::vector <double> y2 = computeYCoordEulerMethod(m, dt, aeroForce, h, y_speed);
	unsigned int end_time_computeY = clock();
	std::cout << end_time_computeY - start_time_computeY << std::endl;
	//std::copy(y.begin(), y.end(), std::ostream_iterator<float>(std::cout,"\n"));
	for (int i = 0; i < y.size(); i++) {
		//	std::cout << y1[i] << "   " << y2[i] << std::endl;
	}
	std::cout << y.size() << ' ' << y2.size() << std::endl;

	double alf = 0;
	int sizeOfcoordinates = y.size();

	std::vector <std::pair<double, double>> time = searchTimeSector(Wx, y, h);

	for (int i = 0; i < time.size(); i++) {
		std::cout << time[i].first << " " << time[i].second << std::endl;
	}
	if (isSingleTrajectory) {
		alf = M_PI / 6;
		std::cout.precision(20);
		std::vector <double> xRK = computeCoordinatesRKMethod(m, dt, Wx, aeroForce, startSpeed*cos(alf), y, time, x_speed);
		std::cout << std::endl;
		std::vector <double> zRK = computeCoordinatesRKMethod(m, dt, Wz, aeroForce, startSpeed*sin(alf), y, time, z_speed);
		x = computeCoordinatesEulerMethod(m, dt, Wx, aeroForce, startSpeed*cos(alf), y, time, x_speed);
		z = computeCoordinatesEulerMethod(m, dt, Wz, aeroForce, startSpeed*sin(alf), y, time, z_speed);
		std::cout << x[x.size() - 1] << " " << z[z.size() - 1] << std::endl;
		std::cout << xRK[xRK.size() - 1] << " " << zRK[zRK.size() - 1] << std::endl;
		std::cout << abs(xRK[xRK.size() - 1] - x[x.size() - 1]) << " " << abs(zRK[zRK.size() - 1] - z[z.size() - 1]) << std::endl;
		writeFile(x, "xTrajectory_Euler.txt");
		writeFile(xRK, "xTrajectory_RK.txt");
	}
	else {
		double  M[89][1];

		for (int i = 0; i < 90; i++) {
			if (i == 89) {
				//	system("pause");
			}
			alf = alf + 2 * M_PI / 360;
			x = computeCoordinatesEulerMethod(m, dt, Wx, aeroForce, startSpeed*cos(alf), y, time, x_speed);
			z = computeCoordinatesEulerMethod(m, dt, Wz, aeroForce, startSpeed*sin(alf), y, time, z_speed);
			x_speed.clear();
			z_speed.clear();
			std::cout << alf << ' ' << startSpeed*cos(alf) << ' ' << x[x.size() - 1] << ' ' << z[z.size() - 1] << ' ' << std::endl;
			replacePath(x, x[x.size() - 1]);
			replacePath(z, z[z.size() - 1]);
			M[i][0] = x[0]; M[i][1] = z[0];
			std::cout << M[i][0] << ' ' << M[i][1] << ' ' << std::endl;
			std::cout << i << std::endl;
		}
	}
	std::ofstream out;
	out.open("\trajectory.txt");
	if (out.is_open())
	{
		for (int i = 0; i < sizeOfcoordinates; i++) {
			//	out << dt*i << ' ' << sqrt(x_speed[i] * x_speed[i] + y_speed[i] * y_speed[i] + z_speed[i] * z_speed[i]) << ' ' << x[i] << ' ' << z[i] << ' ' << y[i] << std::endl;
		}
	}
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	std::cout << search_time << std::endl;

	return 0;
}