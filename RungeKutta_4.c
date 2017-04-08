/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   RungeKutta.c
 * Author: Dominic Gastaldo
 *
 * Created on April 7, 2017, 6:18 PM
 */

#include <stdio.h>
#include <stdlib.h>

/*
 *
 */

double function(double t, double y)
{
    return y;
}

double k1(double t, double y)
{
    return function(t,y);
}

double k2(double t, double y, double h)
{
    return function(t + h/2.0, y + h/2.0 * k1(t,y));
}

double k3(double t, double y, double h)
{
    return function(t + h/2.0, y + h/2.0 * k2(t,y,h));
}

double k4(double t, double y, double h)
{
    return function(t + h, y + h*k3(t,y,h));
}

double RK4(double h, double y, double t, int iterations)
{
    
    
    
    for(int c = 0; c < iterations; c++)
    {
    
        //at nth step
        //RK4 step
        y = y + h/6.0 * (k1(t,y) + 2.0*k2(t,y,h) + 2.0*k3(t,y,h) + k4(t,y,h));
        printf("y = %lf\n", y);
        t = t + h;
        printf("t = %lf\n", t);
    
        //at (n+1)th step
        
    }
}



int main(int argc, char** argv) {
    
    int s = 5;
    
    printf("%d \n",s);
    
    
    RK4(0.001, 1.000, 0.000, 100);
    
    

    return (EXIT_SUCCESS);
}

