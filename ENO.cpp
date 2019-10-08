#include <iostream>
#include <vector>
#include <queue>
#include <math.h>

using namespace std;

/**
 * Inputs
 * curr_points: A vector of the values at the current level in the divided difference table
 * eno_points: A vector containing the given points
 * start_index: The index of the first coordinate in the current divided difference calculation
 * coeffs: A vector containing the coefficients of the interpolating polynomial
 * level: The current level at which to calculate divided differences
 * delta: Used to keep track of which x-coordinates to subtract during divided difference calculations
 */
void divided_difference(queue<float>* curr_points, vector<vector<float> >* eno_points, int start_index, vector<float>* coeffs, int level, int delta)
{
    queue<float>& ref_curr = *curr_points;
    vector<vector<float> >& ref_points = *eno_points;

    coeffs->push_back(ref_curr.front());
    for(int i = 0; i < level - 1; ++i)
    {
        float a, b;
        a = ref_curr.front();
        ref_curr.pop();
        b = ref_curr.front();
        ref_curr.push( (b - a)/(ref_points[start_index + delta][0] - ref_points[start_index][0]) );
        start_index++;
    }
    ref_curr.pop();
}

/**
 * Inputs
 * n: The number of points
 * points: The data points given (x0, y0),...,(x_n, y_n); x0 < ... < x_n
 * x: Some value of x at which to interpolate between x0 and x_n
 * i: The indices x_i and x_i+1 of the points around x
 * p >= 2: The number of nodes to use
 * 
 * Output
 * A vector of the coefficients for the ENO interpolating polynomial
 */
vector<float> eno(int n, vector<vector<float> > points, float x, int i,  int p)
{
    // Do a binary search to find the index of the x-coordinate closest to x
    int mid, left, right;
    left = 0;
    right = n-1;

    cout << "Left: " << left << "  Right: " << right << endl;

    while(left < right)
    {
        mid = (left + right)/2;
        if(points[mid][0] > x)
        {
            right = mid - 1;
        }
        else if(points[mid][0] < x)
        {
            left = mid + 1;
        }
        else
        {
            left = mid;
            right = mid;
        }
        cout << "Left: " << left << "  Mid: " << mid << "  Right: " << right << endl;
    }

    cout << "Found point closest to " << x << ": (" << points[left][0] << ", " << points[left][1] << ")" << endl;
    
    bool is_left;
    // Check if current point found is to the left of x or to the right of x
    if(points[left][0] > x)
    {
        cout << "Point is to the left of " << x << endl;
        is_left = false;
    }
    else
    {
        cout << "Point is to the right of " << x << endl;
        is_left = true;
    }

    // Make sure the left and right-most bounds are within the array
    int far_left, far_right;
    if(is_left)
    {
        far_left = 0 ? (left - (p-2) < 0) : (left - (p-2));
        far_right = n - 1 ? (left + (p-1) < n) : (left + (p-1));
    }
    else
    {
        far_left = 0 ? (left - (p-1) < 0) : (left - (p-1));
        far_right = n - 1 ? (left + (p-2) < n) : (left + (p-2));
    }

    int end = far_left + p;
    vector<float> min_coeffs(p, numeric_limits<float>::max());
    //Loop over all possible groups of p points around x
    while(end != far_right)
    {
        queue<float> curr_points;
        vector<float> coeffs(p);
        int level = p;
        int delta = 1;
        for(int i = far_left; i < end; ++i)
        {
            curr_points.push(points[i][1]);
        }
        while(level > 1)
        {
            divided_difference(&curr_points, &points, far_left, &coeffs, level, delta);
            level--;
            delta++;
        }
        coeffs[p] = curr_points.front(); //The last value in the divided difference table (and the queue)
        if(fabs(coeffs[p]) < fabs(min_coeffs[p])) min_coeffs = coeffs;
        far_left++;
        end++;
    }

    return min_coeffs;
}

int main()
{
    vector<vector<float> > points(4);
    points[0] = vector<float>{-1.0f, 1.0f};
    points[1] = vector<float>{0.0f, 1.0f};
    points[2] = vector<float>{1.0f, 0.0f};
    points[3] = vector<float>{2.0f, 2.0f};

    // vector<float> eno(int n, vector<vector<float> > points, float x, int i,  int p)
    vector<float> coeffs = eno(4, points, 0.5, 1, 2);
    for(auto iter = coeffs.begin(); iter != coeffs.end(); ++iter)
    {
        cout << *iter << "  ";
    }
}