#include "TwoViewTransform.h"

void TVT::FindFund()
{
	// Count Fundamental matrix
	Mat AF = (Mat_<float>(N, 9));
	for (int j = 0; j < N; j++)
	{
		AF.at<float>(j, 0) = points2[j].x * points1[j].x;
		AF.at<float>(j, 1) = points2[j].x * points1[j].y;
		AF.at<float>(j, 2) = points2[j].x;
		AF.at<float>(j, 3) = points2[j].y * points1[j].x;
		AF.at<float>(j, 4) = points2[j].y * points1[j].y;
		AF.at<float>(j, 5) = points2[j].y;
		AF.at<float>(j, 6) = points1[j].x;
		AF.at<float>(j, 7) = points1[j].y;
		AF.at<float>(j, 8) = 1;
	}

	Mat w, u, ut, vt, v;
	SVD::compute(AF, w, u, vt);
	F = (Mat_<float>(3, 3));
	F.at<float>(0, 0) = vt.at<float>(8, 0);
	F.at<float>(0, 1) = vt.at<float>(8, 1);
	F.at<float>(0, 2) = vt.at<float>(8, 2);
	F.at<float>(1, 0) = vt.at<float>(8, 3);
	F.at<float>(1, 1) = vt.at<float>(8, 4);
	F.at<float>(1, 2) = vt.at<float>(8, 5);
	F.at<float>(2, 0) = vt.at<float>(8, 6);
	F.at<float>(2, 1) = vt.at<float>(8, 7);
	F.at<float>(2, 2) = vt.at<float>(8, 8);
}

Mat TVT::FundMatrix()
{
	return F;
}

void TVT::FindProj()
{
	Mat Z = (Mat_<float>(3, 3) << 0, 1, 0, -1, 0, 0, 0, 0, 0);
	Mat W = (Mat_<float>(3, 3) << 0, -1, 0, 1, 0, 0, 0, 0, 1);
	K = (Mat_<float>(3, 3) << f, 0, -wp / 2, 0, f, -hp / 2, 0, 0, 1);
	Mat E = (Mat_<float>(3, 3)); // Essential matrix
	Mat Kt, Wt;
	transpose(K, Kt);
	transpose(W, Wt);
	E = Kt * F * K;
	Mat w, u, ut, vt, v;
	SVD::compute(E, w, u, vt);
	transpose(vt, v);
	transpose(u, ut);
	Mat D = (Mat_<float>(3, 3) << 0, w.at<float>(0), 0, 0, w.at<float>(1), 0, 0, w.at<float>(2), 0);
	Mat R1 = (Mat_<float>(3, 3)); // Rotation matrix
	Mat R2 = (Mat_<float>(3, 3));
	Mat t = (Mat_<float>(1, 3) << u.at<float>(0, 2), u.at<float>(1, 2), u.at<float>(2, 2)); // Translation vector
	for (int i = 0; i < 4; i++)
		P[i] = (Mat_<float>(3, 4));
	R1 = u * W * vt;
	R2 = u * Wt * vt;
	for (int i = 0; i < 3; i++)
	{
		P[0].at<float>(i, 3) = t.at<float>(i);
		P[1].at<float>(i, 3) = -t.at<float>(i);
		P[2].at<float>(i, 3) = t.at<float>(i);
		P[3].at<float>(i, 3) = -t.at<float>(i);
		for (int j = 0; j < 3; j++)
		{
			P[0].at<float>(i, j) = R1.at<float>(i, j);
			P[1].at<float>(i, j) = R1.at<float>(i, j);
			P[2].at<float>(i, j) = R2.at<float>(i, j);
			P[3].at<float>(i, j) = R2.at<float>(i, j);
		}
	}
}

Mat TVT::ProjMatrix()
{
	return P[method];
}

Mat TVT::Compute3Dpts()
{
	Mat P0 = (Mat_<float>(3, 4) << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0);
	Mat X[4];
	for (int i = 0; i < 4; i++)
		X[i] = (Mat_<float>(4, 31));
	for (int i = 0; i < 4; i++)
		P[i] = K * P[i];
	P0 = K * P0;
	// Triangulation
	for (int k = 0; k < 4; k++){
		for (int j = 0; j < 31; j++)
		{
			Mat TriCoef = (Mat_<float>(4, 4));
			TriCoef.row(0) = points1[j].x * P0.row(2) - P0.row(0);
			TriCoef.row(1) = points1[j].y * P0.row(2) - P0.row(1);
			TriCoef.row(2) = points2[j].x * P[k].row(2) - P[k].row(0);
			TriCoef.row(3) = points2[j].y * P[k].row(2) - P[k].row(1);
			Mat w, u, ut, vt, v;
			SVD::compute(TriCoef, w, u, vt);
			transpose(vt, v);
			v.col(3).copyTo(X[k].col(j));
		}

		for (int i = 0; i < 4; i++)
		for (int j = 0; j < 31; j++)
		{
			X[k].at<float>(i, j) = X[k].at<float>(i, j) / X[k].at<float>(3, j);
		}
	}

	// Choose right method
	for (int k = 0; k < 4; k++)
	{
		double mn, mx;
		minMaxIdx(X[k].row(2), &mn, &mx);
		if (mn > 0 && mn * mx > 0)
			method = k;
	}
	return X[method];
}