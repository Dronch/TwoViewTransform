#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/core/core.hpp"

using namespace cv;
using namespace std;

class TVT
{
public:

	TVT(vector<Point2f> p1, vector<Point2f> p2, int focus, int width, int height)
	{
		points1 = p1;
		points2 = p2;
		N = points1.size();
		f = focus;
		wp = width;
		hp = height;
	}

	void FindFund();
	Mat FundMatrix();
	void FindProj();
	Mat ProjMatrix();
	Mat Compute3Dpts();
	int RightMethod();

private:
	vector<Point2f> points1;
	vector<Point2f> points2; // Points from 1st and 2nd pic
	int N; // Number of points
	Mat F; // Fundamental matrix
	Mat P[4]; // Projection matrix
	int method; // choosen method
	Mat K; // Internal camera parameters
	int f, wp, hp; // Camera focus and pixel width and height of picture
};