#include "../inc/geo_util.hpp"

double findTetraAspectRatio(std::vector<std::vector<double>> &_verts)
{
    double in_radius = findInRadius(_verts);

    double circum_radius = findCircumRadius(_verts);

    double l0, l1, l2, l3, l4, l5, max_l;
    l0 = findLength(_verts[0], _verts[1]);
    l1 = findLength(_verts[0], _verts[2]);
    l2 = findLength(_verts[0], _verts[3]);
    l3 = findLength(_verts[1], _verts[2]);
    l4 = findLength(_verts[1], _verts[3]);
    l5 = findLength(_verts[2], _verts[3]);
    max_l = std::max({l0, l1, l2, l3, l4, l5});

    std::cout << in_radius << " " << circum_radius << std::endl; 

    //return max_l / (2*sqrt(6)*in_radius);
    return circum_radius/(3*in_radius);
}

double findCircumRadius(std::vector<std::vector<double>> &_verts)
{
    double circum_radius;
    std::vector<double> formulaN = findformulaN(_verts);
    double formulaD = findformulaD(_verts);

    //std::cout << "formulaN " << formulaN[0] << " " << formulaN[1] << " " << formulaN[2] << std::endl;
    //std::cout << "formulaD " << formulaD << std::endl;

    double norm_formulaN = sqrt(formulaN[0]*formulaN[0] + formulaN[1]*formulaN[1] + formulaN[2]*formulaN[2]);
    circum_radius = norm_formulaN / (2*fabs(formulaD));

    return circum_radius;
}

std::vector<double> findCircumCenter(std::vector<std::vector<double>> &_verts)
{
    std::vector<double> circum_center(3, sizeof(double));

    double formulaD = findformulaD(_verts);
    std::vector<double> formulaN =  findformulaN(_verts);
    double p = 1/(2*fabs(formulaD));

    double cx = _verts[0][0] + p*formulaN[0];
    double cy = _verts[0][1] + p*formulaN[1];
    double cz = _verts[0][2] + p*formulaN[2];

    circum_center = {cx, cy, cz};

    return circum_center;
}

double findInRadius(std::vector<std::vector<double>> &_verts)
{
    double in_radius;

    double x1 = _verts[1][0] - _verts[0][0];
    double y1 = _verts[1][1] - _verts[0][1];
    double z1 = _verts[1][2] - _verts[0][2];

    double x2 = _verts[2][0] - _verts[0][0];
    double y2 = _verts[2][1] - _verts[0][1];
    double z2 = _verts[2][2] - _verts[0][2];

    double x3 = _verts[3][0] - _verts[0][0];
    double y3 = _verts[3][1] - _verts[0][1];
    double z3 = _verts[3][2] - _verts[0][2];

    double x4 = _verts[2][0] - _verts[1][0];
    double y4 = _verts[2][1] - _verts[1][1];
    double z4 = _verts[2][2] - _verts[1][2];

    double x5 = _verts[3][0] - _verts[1][0];
    double y5 = _verts[3][1] - _verts[1][1];
    double z5 = _verts[3][2] - _verts[1][2];


    double a00 = x1*x1 + y1*y1 + z1*z1;
    double a01 = x1*x2 + y1*y2 + z1*z2; 
    double a10 = x2*x1 + y2*y1 + z2*z1; 
    double a11 = x2*x2 + y2*y2 + z2*z2;

    double b00 = x2*x2 + y2*y2 + z2*z2;
    double b01 = x2*x3 + y2*y3 + z2*z3; 
    double b10 = x3*x2 + y3*y2 + z3*z2; 
    double b11 = x3*x3 + y3*y3 + z3*z3;

    double c00 = x3*x3 + y3*y3 + z3*z3;
    double c01 = x3*x1 + y3*y1 + z3*z1; 
    double c10 = x1*x3 + y1*y3 + z1*z3; 
    double c11 = x1*x1 + y1*y1 + z1*z1;

    double d00 = x4*x4 + y4*y4 + z4*z4;
    double d01 = x4*x5 + y4*y5 + z4*z5; 
    double d10 = x5*x4 + y5*y4 + z5*z4; 
    double d11 = x5*x5 + y5*y5 + z5*z5;


    double a1 = sqrt(fabs(a00*a11 - a10*a01))*0.5;
    double a2 = sqrt(fabs(b00*b11 - b10*b01))*0.5;
    double a3 = sqrt(fabs(c00*c11 - c10*c01))*0.5;
    double a4 = sqrt(fabs(d00*d11 - d10*d01))*0.5;

    double a = fabs(a1) + fabs(a2) + fabs(a3) + fabs(a4);
    double vol = findformulaD(_verts)/6;
    in_radius = (3*vol)/a;

    return in_radius;
}

std::vector<double> findInCenter(std::vector<std::vector<double>> &_verts)
{
    std::vector<double> in_center(3, sizeof(double));

    double x1 = _verts[1][0] - _verts[0][0];
    double y1 = _verts[1][1] - _verts[0][1];
    double z1 = _verts[1][2] - _verts[0][2];

    double x2 = _verts[2][0] - _verts[0][0];
    double y2 = _verts[2][1] - _verts[0][1];
    double z2 = _verts[2][2] - _verts[0][2];

    double x3 = _verts[3][0] - _verts[0][0];
    double y3 = _verts[3][1] - _verts[0][1];
    double z3 = _verts[3][2] - _verts[0][2];

    double x4 = _verts[2][0] - _verts[1][0];
    double y4 = _verts[2][1] - _verts[1][1];
    double z4 = _verts[2][2] - _verts[1][2];

    double x5 = _verts[3][0] - _verts[1][0];
    double y5 = _verts[3][1] - _verts[1][1];
    double z5 = _verts[3][2] - _verts[1][2];


    double a00 = x1*x1 + y1*y1 + z1*z1;
    double a01 = x1*x2 + y1*y2 + z1*z2; 
    double a10 = x2*x1 + y2*y1 + z2*z1; 
    double a11 = x2*x2 + y2*y2 + z2*z2;

    double b00 = x2*x2 + y2*y2 + z2*z2;
    double b01 = x2*x3 + y2*y3 + z2*z3; 
    double b10 = x3*x2 + y3*y2 + z3*z2; 
    double b11 = x3*x3 + y3*y3 + z3*z3;

    double c00 = x3*x3 + y3*y3 + z3*z3;
    double c01 = x3*x1 + y3*y1 + z3*z1; 
    double c10 = x1*x3 + y1*y3 + z1*z3; 
    double c11 = x1*x1 + y1*y1 + z1*z1;

    double d00 = x4*x4 + y4*y4 + z4*z4;
    double d01 = x4*x5 + y4*y5 + z4*z5; 
    double d10 = x5*x4 + y5*y4 + z5*z4; 
    double d11 = x5*x5 + y5*y5 + z5*z5;


    double a1 = sqrt(fabs(a00*a11 - a10*a01))*0.5;
    double a2 = sqrt(fabs(b00*b11 - b10*b01))*0.5;
    double a3 = sqrt(fabs(c00*c11 - c10*c01))*0.5;
    double a4 = sqrt(fabs(d00*d11 - d10*d01))*0.5;

    double ix = (a1*_verts[3][0] + a2*_verts[1][0] + a3*_verts[2][0] + a4*_verts[0][0])/(a1+a2+a3+a4);
    double iy = (a1*_verts[3][1] + a2*_verts[1][1] + a3*_verts[2][1] + a4*_verts[0][1])/(a1+a2+a3+a4);
    double iz = (a1*_verts[3][2] + a2*_verts[1][2] + a3*_verts[2][2] + a4*_verts[0][2])/(a1+a2+a3+a4);

    in_center = {ix, iy, iz};

    return in_center;
}

std::vector<double> findformulaN(std::vector<std::vector<double>> &_verts)
{
    std::vector<double> formulaN(3, sizeof(double));

    double x1 = _verts[1][0] - _verts[0][0];
    double y1 = _verts[1][1] - _verts[0][1];
    double z1 = _verts[1][2] - _verts[0][2];

    double x2 = _verts[2][0] - _verts[0][0];
    double y2 = _verts[2][1] - _verts[0][1];
    double z2 = _verts[2][2] - _verts[0][2];

    double x3 = _verts[3][0] - _verts[0][0];
    double y3 = _verts[3][1] - _verts[0][1];
    double z3 = _verts[3][2] - _verts[0][2];

    double a = findLength(_verts[0], _verts[1]);
    double b = findLength(_verts[0], _verts[2]);
    double c = findLength(_verts[0], _verts[3]);

    formulaN[0] = a*a*(y2*z3-z2*y3) + b*b*(y3*z1-z3*y1) + c*c*(y1*z2-z1*y2);
    formulaN[1] = a*a*(z2*x3-x2*z3) + b*b*(z3*x1-x3*z1) + c*c*(z1*x2-x1*z2);
    formulaN[2] = a*a*(x2*y3-y2*x3) + b*b*(x3*y1-y3*x1) + c*c*(x1*y2-y1*x2);

    return formulaN;
}

double findformulaD(std::vector<std::vector<double>> &_verts)
{
    double x1 = _verts[1][0] - _verts[0][0];
    double y1 = _verts[1][1] - _verts[0][1];
    double z1 = _verts[1][2] - _verts[0][2];

    double x2 = _verts[2][0] - _verts[0][0];
    double y2 = _verts[2][1] - _verts[0][1];
    double z2 = _verts[2][2] - _verts[0][2];

    double x3 = _verts[3][0] - _verts[0][0];
    double y3 = _verts[3][1] - _verts[0][1];
    double z3 = _verts[3][2] - _verts[0][2];

    double det = x1*y2*z3 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 - x3*y2*z1;

    //std::cout << x1 << " " <<  y1 << " " << z1 <<std::endl;
    //std::cout << x2 << " " <<  y2 << " " << z2 <<std::endl;
    //std::cout << x3 << " " <<  y3 << " " << z3 <<std::endl;
    //std::cout << "DET" << det << std::endl;

    return det;
}

double findLength(std::vector<double> &_p1, std::vector<double> &_p2)
{
    double length;
    length = sqrt((_p1[0]-_p2[0])*(_p1[0]-_p2[0]) 
                + (_p1[1]-_p2[1])*(_p1[1]-_p2[1]) 
                + (_p1[2]-_p2[2])*(_p1[2]-_p2[2]));
    
    return length;
}

void outputAspectRatio(const std::string &data_dir,
     std::vector<int> &marked_elem, 
     std::vector<std::vector<int>> &cny, 
     std::vector<std::vector<double>> &coor, int nelem_marked)
{
    std::ofstream file("aspectratio.dat");

    std::vector<std::vector<double>> verts(4, std::vector<double>(3));
    for (int i=0; i<nelem_marked; ++i)
    {
        for (int j=0; j<4; ++j)
        {
            for (int k=0; k<3; ++k)
            {
                verts[j][k] = coor[cny[marked_elem[i]][j]][k];
            }
        }
        //std::cout << verts[0][0] << verts[0][1] << verts[0][2] << std::endl;
        //std::cout << verts[1][0] << verts[1][1] << verts[1][2] << std::endl;
        //std::cout << verts[2][0] << verts[2][1] << verts[2][2] << std::endl;
        //std::cout << verts[3][0] << verts[3][1] << verts[3][2] << std::endl;
        //std::cout << std::endl;
        double aspect_ratio = findTetraAspectRatio(verts);
        file << aspect_ratio << std::endl;

    }
    file.close();
}