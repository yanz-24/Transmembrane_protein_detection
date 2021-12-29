import numpy as np
def translate_to_origin():
        return 0
        
def get_rotation_matrix(vec1, vec2):
        """ Find the rotation matrix that aligns vec1 to vec2
        :param vec1: A 3d "source" vector
        :param vec2: A 3d "destination" vector
        :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.

        https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
        """
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix

z_axis = np.array([0, 0, 1])
best_vector = np.array([92.25634277022061, 18.234639485677082, 30.08254356882667])
vec = np.array([1, 0, 0])
mat = get_rotation_matrix(vec1=best_vector, vec2=z_axis)

print(z_axis)
print(np.array([z_axis]).T)
print(mat)
print(np.dot(mat, best_vector))
