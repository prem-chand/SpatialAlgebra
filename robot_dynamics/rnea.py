# from robot_dynamics.rnea import MultiBodySystem, RigidBodyParams
import numpy as np
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class RigidBodyParams:
    """Parameters for a single rigid body/link in the robot"""
    mass: float                    # Mass of the link
    inertia: np.ndarray           # 3x3 inertia tensor
    com_position: np.ndarray      # Center of mass position in link frame
    joint_axis: np.ndarray        # Joint axis in link frame
    parent_to_child: np.ndarray   # Homogeneous transform from parent to child frame


class MultiBodySystem:
    def __init__(self, bodies: List[RigidBodyParams]):
        self.bodies = bodies
        self.num_links = len(bodies)

    def compute_coriolis_terms(self, q: np.ndarray, qdot: np.ndarray) -> np.ndarray:
        """
        Compute Coriolis terms using the Recursive Newton-Euler Algorithm

        Args:
            q: Joint positions (n x 1)
            qdot: Joint velocities (n x 1)

        Returns:
            c: Coriolis force vector (n x 1)
        """
        # Initialize arrays for velocities and accelerations
        v = [np.zeros(3) for _ in range(self.num_links)
             ]        # Linear velocities
        w = [np.zeros(3) for _ in range(self.num_links)
             ]        # Angular velocities
        vdot = [np.zeros(3) for _ in range(self.num_links)
                ]     # Linear accelerations
        wdot = [np.zeros(3) for _ in range(self.num_links)
                ]     # Angular accelerations

        # Forward recursion - Kinematics
        for i in range(self.num_links):
            R = self.bodies[i].parent_to_child[:3,
                                               :3]          # Rotation matrix
            # Position vector
            p = self.bodies[i].parent_to_child[:3, 3]

            if i == 0:
                # Base link
                w[i] = qdot[i] * self.bodies[i].joint_axis
                v[i] = np.zeros(3)
                wdot[i] = np.zeros(3)
                vdot[i] = np.zeros(3)
            else:
                # Propagate velocities from parent to current link
                w[i] = R.T @ w[i-1] + qdot[i] * self.bodies[i].joint_axis
                v[i] = R.T @ (v[i-1] + np.cross(w[i-1], p))

                # Compute accelerations
                wdot[i] = R.T @ wdot[i-1] + \
                    np.cross(R.T @ w[i-1], qdot[i] * self.bodies[i].joint_axis)
                vdot[i] = R.T @ (vdot[i-1] + np.cross(wdot[i-1], p) +
                                 np.cross(w[i-1], np.cross(w[i-1], p)))

        # Backward recursion - Forces and moments
        f = [np.zeros(3) for _ in range(self.num_links)]        # Forces
        n = [np.zeros(3) for _ in range(self.num_links)]        # Moments
        tau = np.zeros(self.num_links)                          # Joint torques

        for i in range(self.num_links - 1, -1, -1):
            com = self.bodies[i].com_position
            m = self.bodies[i].mass
            I = self.bodies[i].inertia

            # Compute force and moment at center of mass
            f[i] = m * (vdot[i] + np.cross(wdot[i], com) +
                        np.cross(w[i], np.cross(w[i], com)))
            n[i] = I @ wdot[i] + np.cross(w[i], I @ w[i])

            if i < self.num_links - 1:
                # Propagate forces and moments from child to parent
                R_next = self.bodies[i+1].parent_to_child[:3, :3]
                p_next = self.bodies[i+1].parent_to_child[:3, 3]

                f[i] += R_next @ f[i+1]
                n[i] += R_next @ n[i+1] + np.cross(p_next, R_next @ f[i+1])

            # Compute joint torque (projection onto joint axis)
            tau[i] = n[i] @ self.bodies[i].joint_axis

        return tau


def skew(v: np.ndarray) -> np.ndarray:
    """Convert vector to skew-symmetric matrix"""
    return np.array([[0, -v[2], v[1]],
                    [v[2], 0, -v[0]],
                    [-v[1], v[0], 0]])


# Example usage

# Create parameters for each link
link1 = RigidBodyParams(
    mass=1.0,
    inertia=np.eye(3),
    com_position=np.array([0, 0, 0.5]),
    joint_axis=np.array([0, 0, 1]),
    parent_to_child=np.eye(4)  # Replace with actual transform
)

# Create the system
robot = MultiBodySystem([link1, ...])  # Add all links

# Compute Coriolis terms
q = np.array([0.1, 0.2, 0.3])      # Joint positions
qdot = np.array([0.1, 0.1, 0.1])   # Joint velocities
coriolis = robot.compute_coriolis_terms(q, qdot)
