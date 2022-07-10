import numpy as np


class DeltaRobot:
    def __init__(self):
        self.sb = 567
        self.sp = 76
        self.L = 524
        self.l = 1244
        self.h = 131
        self.wb = 164
        self.ub = 327
        self.wp = 22
        self.up = 44

        self.B1 = np.array([0, -self.wb, 0])
        self.B2 = np.array([(np.sqrt(3)/2)*self.wb, (1/2)*self.wb, 0])
        self.B3 = np.array([(-np.sqrt(3) / 2) * self.wb, (1 / 2) * self.wb, 0])

    # def __init__(self, R, r, a, b):
    #     self.R = R
    #     self.r = r
    #     self.a = a
    #     self.b = b
    #     self.phi1 = 0
    #     self.phi2 = 2*np.pi/3
    #     self.phi3 = -2 * np.pi / 3

    # def compute_ik_jacobian(self, theta):
    #     th11, th21, th31 = theta[0][0:3]
    #     th12, th22, th32 = theta[1][0:3]
    #     th13, th23, th33 = theta[2][0:3]
    #     j1x = np.sin(th31)*np.cos(th21+th11)*np.cos(self.phi1)+np.cos(th31)*np.sin(self.phi1)
    #     j1y = -np.sin(th31) * np.cos(th21 + th11) * np.sin(self.phi1) + np.cos(th31) * np.cos(self.phi1)
    #     j1z = np.sin(th31)*np.sin(th21+th11)
    #
    #     j2x = np.sin(th32) * np.cos(th22 + th12) * np.cos(self.phi2) + np.cos(th32) * np.sin(self.phi2)
    #     j2y = -np.sin(th32) * np.cos(th22 + th12) * np.sin(self.phi2) + np.cos(th32) * np.cos(self.phi2)
    #     j2z = np.sin(th32) * np.sin(th22 + th12)
    #
    #     j3x = np.sin(th33) * np.cos(th23 + th13) * np.cos(self.phi3) + np.cos(th33) * np.sin(self.phi3)
    #     j3y = -np.sin(th33) * np.cos(th23 + th13) * np.sin(self.phi3) + np.cos(th33) * np.cos(self.phi3)
    #     j3z = np.sin(th33) * np.sin(th23 + th13)
    #
    #     Jp = np.array([[j1x, j1y, j1z], [j2x, j2y, j2z], [j3x, j3y, j3z]])
    #
    #     jth1 = self.a*np.sin(th21)*np.sin(th31)
    #     jth2 = self.a*np.sin(th22) * np.sin(th32)
    #     jth3 = self.a*np.sin(th23) * np.sin(th33)
    #
    #     Jth = np.array([[jth1, 0, 0], [0, jth2, 0], [0, 0, jth3]])
    #
    #     J = np.matmul(np.linalg.inv(Jth), Jp)
    #
    #     return J

    def compute_ik(self, x, y, z):
        a = self.wb-self.up
        b = self.sp/2-np.sqrt(3)*self.wb/2
        c = self.wp-self.wb/2
        E1 = 2 * self.L * (y + a)
        F1 = 2 * z * self.L
        G1 = x ** 2 + y ** 2 + z ** 2 + a ** 2 +self.L ** 2 + 2 * y * a - self.l ** 2
        E2 = -self.L*(np.sqrt(3)*(x+b)+y+c)
        F2 = 2*z*self.L
        G2 = x ** 2 + y ** 2 + z ** 2 + b**2 + c**2 + self.L**2 + 2*(x*b+y*c)-self.l**2
        E3 = self.L*(np.sqrt(3)*(x-b)-y-c)
        F3 = 2*z*self.L
        G3 = x ** 2 + y ** 2 + z ** 2 + b**2 + c**2 + self.L**2 + 2*(-x*b+y*c)-self.l**2
        sign = 1
        t1 = (-F1 + sign * np.sqrt(E1 ** 2 + F1 ** 2 - G1 ** 2)) / (G1 - E1)
        t2 = (-F2 + sign * np.sqrt(E2 ** 2 + F2 ** 2 - G2 ** 2)) / (G2 - E2)
        t3 = (-F3 + sign * np.sqrt(E3 ** 2 + F3 ** 2 - G3 ** 2)) / (G3 - E3)
        theta_1 = 2 * np.arctan(t1)
        theta_2 = 2 * np.arctan(t2)
        theta_3 = 2 * np.arctan(t3)
        theta = np.array([theta_1, theta_2, theta_3])

        # val1 = (G1-E1)*t1**2+(2*F1)*t1+(G1+E1)
        # val2 = (G2 - E2) * t2 ** 2 + (2 * F2) * t2 + (G2 + E2)
        # val3 = (G3 - E3) * t3 ** 2 + (2 * F3) * t3 + (G3 + E3)
        # print('here1', val1, val2, val3)
        #
        # val4 = E1*np.cos(theta_1)+F1*np.sin(theta_1)+G1
        # val5 = E2 * np.cos(theta_2) + F2 * np.sin(theta_2) + G2
        # val6 = E3 * np.cos(theta_3) + F3 * np.sin(theta_3) + G3
        # print('here2', val4, val5, val6)
        #
        # val9 = self.L*(np.sqrt(3)*(x-b)-y-c) * np.cos(theta_3) + 2*z*self.L * np.sin(theta_3) + x ** 2 + y ** 2 + z ** 2 + b**2 + c**2 + self.L**2 + 2*(-x*b+y*c)-self.l**2
        # print('xyz', x, y, z)
        # print('here3', val9)


        return theta

    def forward_ik(self, theta):
        theta_1, theta_2, theta_3 = theta

        A1v = np.array([0, -self.wb-self.L*np.cos(theta_1)+self.up, -self.L*np.sin(theta_1)])
        A2v = np.array([(np.sqrt(3)/2)*(self.wb+self.L*np.cos(theta_2))-self.sp/2, (1/2)*(self.wb+self.L*np.cos(theta_2))-self.wp, -self.L*np.sin(theta_2)])
        A3v = np.array([(-np.sqrt(3)/2)*(self.wb+self.L*np.cos(theta_3))+self.sp/2, (1/2)*(self.wb+self.L*np.cos(theta_3))-self.wp, -self.L*np.sin(theta_3)])

        x, y, z = self.trilaterate(A1v, A2v, A3v)

        return x, y, z

    def sphere_intersection(self, A1, A2, A3):
        x1, y1, z1 = A1
        x2, y2, z2 = A2
        x3, y3, z3 = A3
        zn = z1
        r1, r2, r3 = 3*[self.l]

        a = 2 * (x3-x1)
        b = 2*(y3-y1)
        c = r1**2-r3**2-x1**2-y1**2+x3**2+y3**2
        d = 2*(x3-x2)
        e = 2*(y3-y2)
        f = r2**2-r3**2-x2**2-y2**2+x3**2+y3**2

        x = (c*e-b*f)/(a*e-b*d)
        y = (a*f - c*d) / (a*e - b * d)

        A = 1
        B = -2*zn
        C = zn**2-r1**2+(x-x1)**2+(y-y1)**2

        zp = (-B + np.sqrt(B ** 2 - (4 * A * C))) / (2 * A)
        zm = (-B - np.sqrt(B ** 2 - (4 * A * C))) / (2 * A)
        z = min(zp, zm)

        val1 = (x-x1)**2+(y-y1)**2+(z-zn)**2-r1**2
        val2 = (x - x2) ** 2 + (y - y2) ** 2 + (z - zn) ** 2 - r2 ** 2
        val3 = (x - x3) ** 2 + (y - y3) ** 2 + (z - zn) ** 2 - r3 ** 2

        print('here', val1, val2, val3)

        return x, y, z

    def sphere_intersection2(self, A1, A2, A3):
        x1, y1, z1 = A1
        x2, y2, z2 = A2
        x3, y3, z3 = A3
        r1, r2, r3 = 3*[self.l]

        a11 = 2 * (x3 - x1)
        a12 = 2 * (y3 - y1)
        a13 = 2 * (z3 - z1)
        a21 = 2 * (x3 - x2)
        a22 = 2 * (y3 - y2)
        a23 = 2 * (z3 - z2)
        b1 = r1**2-r3**2-x1**2-y1**2-z1**2+x3**2+y3**2+z3**2
        b2 = r2 ** 2 - r3 ** 2 - x2 ** 2 - y2 ** 2 - z2 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2

        a1 = a11/a13-a21/a23
        a2 = a12/a13-a22/a23
        a3 = b2/a23-b1/a13
        a4 = -a2/a1
        a5 = -a3/a1
        a6 = (-a21*a4-a22)/a23
        a7 = (b2-a21*a5)/a23

        a = a4**2+1+a6**2
        b = 2*a4*(a5-x1)-2*y1+2*a6*(a7-z1)
        c = a5*(a5-2*x1)+a7*(a7-2*z1)+x1**2+y1**2+z1**2-r1**2

        yp = (-b+np.sqrt(b**2-4*a*c))/(2*a)
        ym = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
        xp = a4*yp+a5
        xm = a4 * ym + a5
        zp = a6 * yp + a7
        zm = a6 * ym + a7
        if zm < zp:
            return xm, ym, zm
        else:
            return xp, yp, zp

    def trilaterate(self, A1, A2, A3):
        r1, r2, r3 = 3 * [self.l]
        temp1 = A2 - A1
        e_x = temp1 / np.linalg.norm(temp1)
        temp2 = A3 - A1
        i = np.dot(e_x, temp2)
        temp3 = temp2 - i * e_x
        e_y = temp3 / np.linalg.norm(temp3)
        e_z = np.cross(e_x, e_y)
        d = np.linalg.norm(A2 - A1)
        j = np.dot(e_y, temp2)
        x = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
        y = (r1 * r1 - r3 * r3 - 2 * i * x + i * i + j * j) / (2 * j)
        temp4 = r1 * r1 - x * x - y * y
        if temp4 < 0:
            raise Exception("The three spheres do not intersect!");
        z = np.sqrt(temp4)
        p_12_a = A1 + x * e_x + y * e_y + z * e_z
        p_12_b = A1 + x * e_x + y * e_y - z * e_z
        if p_12_a[2]<p_12_b[2]:
            return p_12_a
        else:
            return p_12_b

    def compute_medial_joints(self, theta):
        theta_1, theta_2, theta_3 = theta
        A1 = np.array([0, -self.wb-self.L*np.cos(theta_1), -self.L*np.sin(theta_1)])
        A2 = np.array([(np.sqrt(3)/2)*(self.wb+self.L*np.cos(theta_2)), (1/2)*(self.wb+self.L*np.cos(theta_2)), -self.L * np.sin(theta_2)])
        A3 = np.array(
            [(-np.sqrt(3) / 2) * (self.wb + self.L * np.cos(theta_3)), (1 / 2) * (self.wb + self.L * np.cos(theta_3)),
             -self.L * np.sin(theta_3)])

        return A1, A2, A3

    def compute_distal_joints(self, x, y, z):
        P = np.array([x, y, z])
        P1 = P + np.array([0, -self.up, 0])
        P2 = P + np.array([self.sp/2, self.wp, 0])
        P3 = P + np.array([-self.sp/2, self.wp, 0])

        return P1, P2, P3

    def distal_link_vector(self, theta, x, y, z):
        theta_1, theta_2, theta_3 = theta
        a = self.wb-self.up
        b = (self.sp/2)-(np.sqrt(3)/2)*self.wb
        c = self.wp-(1/2)*self.wb
        l1 = np.array([x, y + self.L*np.cos(theta_1)+a, z + self.L*np.sin(theta_1)])
        l2 = np.array([x-(np.sqrt(3)/2)*self.L*np.cos(theta_2)+b, y -(1/2)*self.L * np.cos(theta_2) + c, z + self.L * np.sin(theta_2)])
        l3 = np.array([x + (np.sqrt(3) / 2) * self.L * np.cos(theta_3) - b, y - (1 / 2) * self.L * np.cos(theta_3) + c,
                       z + self.L * np.sin(theta_3)])

        return l1, l2, l3

    def constraint_equations(self, theta, x, y, z):
        theta_1, theta_2, theta_3 = theta
        a = self.wb - self.up
        b = (self.sp / 2) - (np.sqrt(3) / 2) * self.wb
        c = self.wp - (1 / 2) * self.wb
        eq1 = 2*self.L*(y+a)*np.cos(theta_1)+2*z*self.L*np.sin(theta_1)+x**2+y**2+z**2+a**2+self.L**2+2*y*a-self.l**2
        eq2 = -self.L*(np.sqrt(3)*(x+b)+y+c)*np.cos(theta_2)+2*z*self.L*np.sin(theta_2)+x**2+y**2+z**2+b**2+c**2+self.L**2+2*x*b+2*y*c-self.l**2
        eq3 = self.L * (np.sqrt(3) * (x - b) - y - c) * np.cos(theta_3) + 2 * z * self.L * np.sin(
            theta_3) + x ** 2 + y ** 2 + z ** 2 + b ** 2 + c ** 2 + self.L ** 2 - 2 * x * b + 2 * y * c - self.l ** 2

        # print('xyz', x, y, z)
        eq3b = self.L*(np.sqrt(3)*(x-b)-y-c) * np.cos(theta_3) + 2*z*self.L * np.sin(theta_3) + x ** 2 + y ** 2 + z ** 2 + b**2 + c**2 + self.L**2 + 2*(-x*b+y*c)-self.l**2

        return eq1, eq2, eq3, eq3b

    def compute_jacobian(self, theta):
        x, y, z = self.forward_ik(theta)
        theta_1, theta_2, theta_3 = theta
        a = self.wb-self.up
        b = self.sp/2-np.sqrt(3)*self.wb/2
        c = self.wp-self.wb/2

        # A * dx = B * dtheta (A and B Jacobian matrices here!)
        A1 = [x, y+self.L*np.cos(theta_1), z+self.L*np.sin(theta_1)]
        A2 = [2*(x+b)-np.sqrt(3)*self.L*np.cos(theta_2), 2*(y+c) - self.L * np.cos(theta_2), 2*(z + self.L * np.sin(theta_2))]
        A3 = [2 * (x - b) + np.sqrt(3) * self.L * np.cos(theta_3), 2 * (y + c) - self.L * np.cos(theta_3),
              2 * (z + self.L * np.sin(theta_3))]
        A = np.array([A1, A2, A3])

        b11 = self.L*((y+a)*np.sin(theta_1)-z*np.cos(theta_1))
        b22 = -self.L*(np.sqrt(3)*((x+b)+y+c)*np.sin(theta_2)+2*z*np.cos(theta_2))
        b33 = self.L*(np.sqrt(3)*((x-b)-y-c)*np.sin(theta_3)-2*z*np.cos(theta_3))
        B = np.diag([b11, b22, b33])

        J = np.matmul(np.linalg.inv(A), B)

        return J
