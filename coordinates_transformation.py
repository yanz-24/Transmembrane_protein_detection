import numpy as np

def translate_df(df, vector):
	"""Translate the coordinate system of pdb dataframe with a vector.

	:param df: Pdb pandas dataframe containing 3D coordinates of atoms
	:param vector: translation vector
	:return df: df with translated coordinates
	"""
	df["x_coord"] = df["x_coord"] - vector[0]
	df["y_coord"] = df["y_coord"] - vector[1]
	df["z_coord"] = df["z_coord"] - vector[2]
	return df


def rotate_df(df, matrix):
	"""Rotate the coordinate system of pdb dataframe with a matrix.

	:param df: Pdb pandas dataframe containing 3D coordinates of atoms
	:param vector: rotation matrix (3x3)
	:return df_rotated: df with rotated coordinates
	"""
	df_rotated = df
	df_rotated["x_coord"] = df["x_coord"] * matrix[0,0] + df["y_coord"] * matrix[0,1] + df["z_coord"] * matrix[0,2]
	df_rotated["y_coord"] = df["x_coord"] * matrix[1,0] + df["y_coord"] * matrix[1,1] + df["z_coord"] * matrix[1,2]
	df_rotated["z_coord"] = df["x_coord"] * matrix[2,0] + df["y_coord"] * matrix[2,1] + df["z_coord"] * matrix[2,2]
	return 	df_rotated


def get_rotation_matrix(vec1, vec2):
	""" Find the rotation matrix that aligns vec1 to vec2.

	:param vec1: A 3d "source" vector
	:param vec2: A 3d "destination" vector
	:return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.

	source: https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
	"""
	a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
	v = np.cross(a, b)
	c = np.dot(a, b)
	s = np.linalg.norm(v)
	kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
	rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
	return rotation_matrix
