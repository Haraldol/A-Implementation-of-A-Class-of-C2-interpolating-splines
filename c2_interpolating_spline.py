import pygame
from pygame.locals import *
from sys import exit
import numpy as np
#
# Author: Harald Olin
#
# This program implements two functions for interpolation, Bezier interpolation and circle interpolation
# and utilises the blending function mentioned in Cem Yuksels (2020) to achieve C^2 interpolating splines.   
# It also provides a user iterface where one can place points and choose  which interpolation to use between the points. 
# 
# Cem Yuksel. 2020.
# A Class of C2 Interpolating Splines.
# ACM Trans. Graph. 39, 5, Article 160 (October 2020), 14 pages.
# https://doi.org/10.1145/3400301
#

pygame.init()
screen = pygame.display.set_mode((900, 600))

pygame.display.set_caption("C2 Interpolating spline")

# Initiate Colors
text_color = (255,255,255)
hovering_color = (170,170,170)
button_color = (100,100,100)
menu_color =(128,128,128) 

# Defining a font
smallfont = pygame.font.SysFont('Corbel',30)
# Text for buttons 
bezier_interpolation_text_on = smallfont.render('Bezier Interpolation On', True, text_color)
bezier_interpolation_text_off = smallfont.render('Bezier Interpolation Off', True, text_color)
bezier_blending_text_on = smallfont.render('Bezier Blending On', True, text_color)
bezier_blending_text_off = smallfont.render('Bezier Blending Off', True, text_color)
circle_interpolation_text_on = smallfont.render('Circle Interpolation On', True, text_color)
circle_interpolation_text_off = smallfont.render('Circle Interpolation Off', True, text_color)
circle_blending_text_on = smallfont.render('Circle Blending On', True, text_color)
circle_blending_text_off = smallfont.render('Circle Blending Off', True, text_color)
remove_text = smallfont.render('Remove', True, text_color)

# Global variables 
run = True
point_i = 0
ctrl_points = []
lower_bezier_ctrl_points = [(0,0),(0,0),(0,0)]
upper_bezier_ctrl_points = [(0,0),(0,0),(0,0)]
draw_bezier = False
draw_blend = False
draw_circle = False
draw_blending_circle = False

def quad_bezier(t, control_points):
    """Calculates and returns the quadratic bezier interpolation points"""
    return ( control_points[0][0]*(1-t)**2 + 2*(1-t)*t*control_points[1][0] + control_points[2][0]*t**2,
             control_points[0][1]*(1-t)**2 + 2*(1-t)*t*control_points[1][1] + control_points[2][1]*t**2)

def bezier_blending(theta,lower_bezier_ctrl_points,upper_bezier_ctrl_points,t_i,t_ip1):
    """
    The blending function c_i, using bezier interpolation for the function F_i as mentioned in the paper.
    Following the mathematical formula:
    C(theta) = cos(theta)^2 * F_i(theta + pi/2) + sin(theta)^2 * F_{i+1}(theta)
    where F_i is our chosen interpolation function, in this case, the bezier interpolation.
    The conversion from F_i(theta) to bezier is the following:      
    F_i(theta) = (2*t_i*theta)/pi for  0 <= theta =< pi/2
    F_i(theta) = t_i + (2*(1-t_i)*(theta - pi/2))/pi
    This makes sure that: 
    F_i(0)=B(0)=p_i-1,
    F_i(pi/2)=B(t_i)=p_i,
    F_i(pi)=B(1)=p_i+1,
    """
    theta_i = theta + np.pi/2
    theta_ip1 = theta

    if(np.pi/2 < theta_i):
        theta_i = t_i + (2*(1-t_i)*(theta_i-np.pi/2))/np.pi
    else: 
        theta_i = (2*t_i*theta_i)/np.pi

    if( theta_ip1 <= np.pi/2):
        theta_ip1 = (2*t_ip1*theta_ip1)/np.pi
    else: 
        theta_ip1 = t_ip1 + (2*(1-t_ip1)*(theta_ip1-np.pi/2))/np.pi
        
    f_i  = quad_bezier(theta_i, lower_bezier_ctrl_points)
    f_ip1 = quad_bezier(theta_ip1, upper_bezier_ctrl_points)

    fi = np.array([f_i[0], f_i[1]])
    fip1 = np.array([f_ip1[0], f_ip1[1]])

    c_i = np.cos(theta)**2 * fi + np.sin(theta)**2 * fip1
    return c_i

def find_ti(p_0,p_1,p_2):
    """
    Finds the best suited t_i involved in finding the 2nd
    control point in the bezier curve, 
    by solving the cubic polynom: 
    ||p_{i+1} - p_{i-1} ||^2*t_i^3 + 
    3*(p_{i+1} - p_{i-1}) * (p_{i-1} - p_i)*t_i^2 +
    (3p_{i-1} - 2p_i - p_{i+1}) * (p_{i-1} - p_i)t_i -
    ||p_i-1 - p_i ||^2 = 0
    as proposed in the paper.
    """
    # Transform tuple to np.arrays for mathematical operations
    p0 = np.array([p_0[0], p_0[1]])
    p1 = np.array([p_1[0], p_1[1]])
    p2 = np.array([p_2[0], p_2[1]])
    # Calculate the coefficients a, b, c, d
    a = np.linalg.norm(p2 - p0)**2
    b = 3 * np.dot(p2 - p0, p0 - p1)
    c = np.dot(3 * p0 - 2 * p1 - p2, p0 - p1)
    d = -1*np.linalg.norm(p0 - p1)**2

    coefficients = [a,b,c,d]
    all_roots = np.roots(coefficients)
    # find the one root we are looking for 
    for root in all_roots:
        if np.isreal(root) and 0 <= root.real <= 1:
            t_i = root.real
            return t_i
    
def find_b1(t_i,p_0,p_1,p_2):
    """
    finds the midle control point b_{i,1} for the quadratic bezier curve by solving: 
    b_{i,1} = p_i - (1 - t_i)^2*b_{i,0} - t_i^2 * b_{i,2} / ( 2 * (1 - t_i) * t_i
    as stated in the paper.
    Note that p0 = b_{i,0}, p2 = b_{i,2}
    """
    # Transform tuple to np.arrays for mathematical operations
    p0 = np.array([p_0[0], p_0[1]])
    p1 = np.array([p_1[0], p_1[1]])
    p2 = np.array([p_2[0], p_2[1]])

    numerator = p1-(1-t_i)**2 * p0 - t_i**2 * p2
    denominator = 2 * (1-t_i) * t_i
    b1 = numerator/denominator
    return b1 

def find_circle(p1,p2,p3):
    """
    function to find the circle given 3 points by solving 
    (x_i-x0)^2 + (y_i-y0)^2 -r^2 = 0, where x0,y0 and r are unknown.
    """
    x1,y1 = p1
    x2,y2 = p2
    x3,y3 = p3
    A = np.array ([[2*(x2-x1), 2*(y2-y1)],
                   [2*(x3-x1), 2*(y3-y1)] ])
    
    B = np.array([[x2**2 - x1**2 + y2**2 - y1**2],
                  [x3**2 - x1**2 + y3**2 - y1**2] ])
    
    center = np.linalg.solve(A,B)
    # The coordinates of the center of the circle
    x0,y0 = center.flatten()
    r = np.sqrt((x1-x0)**2 + (y1-y0)**2) 
    return x0,y0,r


def circle_constants(p1,p2,p3,x0,y0):
    """
    calculates the constants that are needed in the circle interpolation,
    using the equation F'(t)= cos(alpha*t+delta)*u + sin(alpha*t+delta)v+q.
    Since we want to solve f(-delta/alpha)=p_i,
    we get F'(-delta/alpha)=1*u+q=p_i --> u = p_i-q,
    But in the paper, this is how v is defined. Therefore 
    I've changed, u = p_i - q, and v to be perpendicular to u.
    In other words, the paper might have the wrong equation.
    """
    center = np.array([x0,y0])
    # Find vector v_i
    u_i = np.array(p2) - center
    #u_i perpendicular to v_i
    v_i = np.array([-u_i[1],u_i[0]])
    #find delta, solve F(0)=p1, and use arctan2 by converting to v_i,u_i coordinates
    delta_i = np.arctan2(np.array(p1 - center).dot(v_i), np.array(p1 - center).dot(u_i))
    #find alpha + delta, solve F(1)=p3, same as above
    alpha_plus_delta = np.arctan2(np.array(p3 - center).dot(v_i), np.array(p3 - center).dot(u_i))
    alpha_i = alpha_plus_delta - delta_i
    return v_i,u_i,alpha_i,delta_i,center

def circle_interpolation(v_i,u_i,alpha,delta,center,t):
    """
    Computes the circle interpolation function:
    F'(t) = cos(alpha * t+delta) * v + sin(alpha*t+delta)*u + q
    """
    theta = alpha * t + delta
    point = center + np.cos(theta) * u_i + np.sin(theta) * v_i
    return point

def circle_blending(v_i,u_i,alpha_i,delta_i,center_i,t_i,
                     v_ip1,u_ip1,alpha_ip1,delta_ip1,center_ip1,t_ip1,theta):
    """
    Run the circle blending function with the circle interpolation function F'.
    Makes the global parametrization F(0)=F'(0)=p_i-1, F(pi/2)=F'(ti)=p_i, F(pi)=F'(1)=p_i+1,
    and calculates the blending function
    C_i = cos(theta)^2*F_i(theta + pi/2) + sin(theta)^2 * F_i+1(theta)
    """
    theta_i = theta + np.pi/2
    theta_ip1 = theta
    if(np.pi/2 < theta_i):
        theta_i = t_i + (2*(1-t_i)*(theta_i-np.pi/2))/np.pi
    else: 
        theta_i = (2*t_i*theta_i)/np.pi

    if( theta_ip1 <= np.pi/2):
        theta_ip1 = (2*t_ip1*theta_ip1)/np.pi
    else: 
        theta_ip1 = t_ip1 + (2*(1-t_ip1)*(theta_ip1-np.pi/2))/np.pi
    
    f_i  = circle_interpolation(v_i,u_i,alpha_i,delta_i,center_i,theta_i)
    f_ip1 = circle_interpolation(v_ip1,u_ip1,alpha_ip1,delta_ip1,center_ip1,theta_ip1)
    fi = np.array([f_i[0], f_i[1]])
    fip1 = np.array([f_ip1[0], f_ip1[1]])
    c_i = np.cos(theta)**2 * fi + np.sin(theta)**2 * fip1
    return c_i


# The main loop for the program
while run:
    for event in pygame.event.get():
        if event.type == QUIT:
            run = False

        # Get the mouse position
        mouse = pygame.mouse.get_pos()

        # if the mouse is in the paintable area
        if 0 <= mouse[0] <= 600 and 0 <= mouse[1]<= 600:
            if event.type == MOUSEBUTTONDOWN:
                ctrl_points.append( pygame.mouse.get_pos())
                pygame.draw.circle(screen, (255, 255, 255), (ctrl_points[point_i][0], ctrl_points[point_i][1]) , 5)
                point_i = (point_i + 1) 

        # Draws Bezier interpolation function
        if 3 <= len(ctrl_points) and draw_bezier:
            n = len(ctrl_points)-1
            t_i = find_ti(ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n])
            b_i_1 = find_b1(t_i,ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n])
            upper_bezier_ctrl_points = [ctrl_points[n-2],b_i_1,ctrl_points[n]]
            lastPoint = None
            # Draw the bezier curves from the 3 latest points.
            lower_bound = 0
            #Interpolate from t_i instead of 0 for a continous line.
            if(3 < len(ctrl_points)):
                lower_bound = t_i
            for t in np.arange(lower_bound, 1, 0.01):
                point = quad_bezier(t,upper_bezier_ctrl_points)
                if lastPoint is not None:
                    pygame.draw.line(screen, (255, 255, 0), lastPoint, point, 1)
                lastPoint = point

        # Draws the Bezier Blending function
        if 4 <= len(ctrl_points) and draw_blend:
            n = len(ctrl_points)-1
            # calculate t_i+1 and b_i+1_1
            t_ip1 = find_ti(ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n])
            b_ip1_1 = find_b1(t_ip1,ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n])
            upper_bezier_ctrl_points = [ctrl_points[n-2],b_ip1_1,ctrl_points[n]]
            # calculate t_i and b_i_1
            t_i = find_ti(ctrl_points[n-3],ctrl_points[n-2],ctrl_points[n-1])
            b_i_1 = find_b1(t_i,ctrl_points[n-3],ctrl_points[n-2],ctrl_points[n-1])
            lower_bezier_ctrl_points = [ctrl_points[n-3],b_i_1,ctrl_points[n-1]]
            lastPoint = None

            for theta in np.arange(0, np.pi/2, np.pi/200):
                point = bezier_blending(theta,lower_bezier_ctrl_points,upper_bezier_ctrl_points,t_i,t_ip1)
                if lastPoint is not None:
                    pygame.draw.line(screen, (255, 0, 0), lastPoint, point, 1)
                lastPoint = point

        # Draws the Circle interpolation function
        if 3 <= len(ctrl_points) and draw_circle:
            n = len(ctrl_points) - 1
            lower_bound = 0
            points =[]
            x0,y0,r = find_circle(ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n])
            v_i,u_i,alpha,delta,center = circle_constants(ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n],x0,y0)
            #Interpolate from t_i instead of 0 for a continous line.
            if 4<=len(ctrl_points):
                lower_bound = -1*delta/alpha
            for t in np.arange(lower_bound, 1.01, 0.01):
                points.append(circle_interpolation(v_i,u_i,alpha,delta,center,t))
            pygame.draw.lines(screen, (0, 255, 0), False, points, 2)
        
        # Draws the Circle blending function
        if 4 <= len(ctrl_points) and draw_blending_circle:
            n = len(ctrl_points) - 1
            points =[]
            # finding the F_i_plus_1 variables
            x0,y0,r = find_circle(ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n])
            v_ip1,u_ip1,alpha_ip1,delta_ip1,center_ip1 = circle_constants(ctrl_points[n-2],ctrl_points[n-1],ctrl_points[n],x0,y0)
            t_ip1 = -1*delta_ip1/alpha_ip1
            # finding the F_i variables
            x0,y0,r = find_circle(ctrl_points[n-3],ctrl_points[n-2],ctrl_points[n-1])
            v_i,u_i,alpha_i,delta_i,center_i = circle_constants(ctrl_points[n-3],ctrl_points[n-2],ctrl_points[n-1],x0,y0)
            t_i = -1*delta_i/alpha_i
            for theta in np.arange(0, np.pi/2,np.pi/200 ):
                points.append(circle_blending(v_i,u_i,alpha_i,delta_i,center_i,t_i,v_ip1,u_ip1,alpha_ip1,delta_ip1,center_ip1,t_ip1,theta))
            pygame.draw.lines(screen, (0, 0, 255), False, points, 2)

        #Background color for the menu
        pygame.draw.rect(screen,menu_color,[600,0,300,600])

        # Bezier interpolating button
        if 625 <= mouse[0] <= 625+250 and 50 <= mouse[1]<= 50+50:
            pygame.draw.rect(screen,hovering_color,[625,50,250,50])
            if event.type == MOUSEBUTTONDOWN:
                if not draw_bezier:
                    draw_bezier = True
                else:
                    draw_bezier = False
        else: 
            pygame.draw.rect(screen,button_color,[625,50,250,50])
        
        # Bezier Blending button
        if 625 <= mouse[0] <= 625+250 and 150 <= mouse[1]<= 150+50:
            pygame.draw.rect(screen,hovering_color,[625,150,250,50])
            if event.type == MOUSEBUTTONDOWN:
                if not draw_blend:
                    draw_blend = True
                else:
                    draw_blend = False
        else: 
            pygame.draw.rect(screen,button_color,[625,150,250,50])

        # Circle interpolation button
        if 625 <= mouse[0] <= 625+250 and 250 <= mouse[1]<= 250+50:
            pygame.draw.rect(screen,hovering_color,[625,250,250,50])
            if event.type == MOUSEBUTTONDOWN:
                if not draw_circle:
                    draw_circle = True
                else:
                    draw_circle = False       
        else:
            pygame.draw.rect(screen,button_color,[625,250,250,50])

        # Circle blending button
        if 625 <= mouse[0] <= 625+250 and 350 <= mouse[1]<= 350+50:
            pygame.draw.rect(screen,hovering_color,[625,350,250,50])
            if event.type == MOUSEBUTTONDOWN:
                if not draw_blending_circle:
                    draw_blending_circle = True
                else:
                    draw_blending_circle = False
        else:
            pygame.draw.rect(screen,button_color,[625,350,250,50])

        # Remove button
        if 625 <= mouse[0] <= 625+250 and 450 <= mouse[1]<= 450+50:
            pygame.draw.rect(screen,hovering_color,[625,450,250,50])
            if event.type == MOUSEBUTTONDOWN:
                screen.fill(0)
                ctrl_points.clear() 
                point_i = 0
        else:
            pygame.draw.rect(screen,button_color,[625,450,250,50])
        
        # Display texts on screen
        # Bezier interpolation text
        if draw_bezier:
            screen.blit(bezier_interpolation_text_on,(580+50,60))    
        else:
            screen.blit(bezier_interpolation_text_off,(580+50,60))
        #Bezier blending text
        if draw_blend:
            screen.blit(bezier_blending_text_on,(580+50,160))
        else:
            screen.blit(bezier_blending_text_off,(580+50,160))
        # Circle interpolation text
        if draw_circle:
            screen.blit(circle_interpolation_text_on,(580+50,260))
        else:
            screen.blit(circle_interpolation_text_off,(580+50,260))
        # Circle blending text
        if draw_blending_circle:
            screen.blit(circle_blending_text_on,(580+50,360))
        else:
            screen.blit(circle_blending_text_off,(580+50,360))
        # Remove text
        screen.blit(remove_text,(650+50,460))

        pygame.display.update()

pygame.quit()
exit()