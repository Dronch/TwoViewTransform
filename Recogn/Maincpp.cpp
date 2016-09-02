#include "TwoViewTransform.h"
#include <iostream>
#include <fstream>


int main(int argc, char* argv[])
{
	ifstream input("input.txt");
	vector<Point2f> points1, points2;
	string line;
	while (getline(input, line))
	{
		stringstream ss(line);
		float x1, y1, x2, y2;
		ss >> x1 >> y1 >> x2 >> y2;
		points1.push_back(Point2f(x1, y1));
		points2.push_back(Point2f(x2, y2));
	}
	TVT *test = new TVT(points1, points2, 3800, 4928, 3264);
	test->FindFund();
	test->FindProj();
	ofstream out3d("3d.txt");
	out3d << test->Compute3Dpts();
	out3d.close();
	ofstream output("output.txt");
	output << "Fundamental matrix" << endl << test->FundMatrix() << endl << endl;
	output << "Projection matrix" << endl << test->ProjMatrix() << endl << endl;
	output << "You can find 3D points in 3d.txt";
	output.close();
	delete test;
	return 0;
}