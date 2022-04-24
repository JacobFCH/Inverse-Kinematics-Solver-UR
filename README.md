# Inverse-Kinematics-Solver-UR

Inverse Kinematics solver for UR Robots, based on the paper written by Rasmus Skovgaard Andersen http://rasmusan.blog.aau.dk/files/ur5_kinematics.pdf

This inverse kinematics solver should work for all UR robots as they have the kinematic structure.

Dependencies can be found in the requirements.txt file and can be installed using the command bellow:
```
pip install -r requirements.txt
```

It is also recomended to use ```np.set_printoptions(suppress=True)``` in order to supress numpy using scientific notation

The IK Solver can be initialized as:

```
a = np.array([0, -0.425, -0.39225, 0, 0, 0])
d = np.array([0.089159, 0, 0, 0.10915, 0.09465, 0.0823])
alpha = np.array([0, np.deg2rad(90), 0, 0, np.deg2rad(90), np.deg2rad(-90)])

ik = ikSolver(a, d, alpha)
```

The DH parameters can be found on the UR website https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/

**Note that alpha is offset by one element this is due to the implementation explained in the paper**

For simplicity the transformation matrix can be created using

```
T = ik.create_Transformation_Matrix(position,orientation)
```

Lastly the IK solver be used to compute the inverse kinematics solutions to the given TCP frame

```
qs = ik.solveIK(T)
```

Returning every inverse kinematics solution

```
q, qs = ik.solveIK(T, last_q)
```

If used with the previous configuration ```last_q``` then the method will return the IK solution that is closest as well as the entire list of IK solutions