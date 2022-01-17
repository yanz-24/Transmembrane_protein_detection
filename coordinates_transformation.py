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
	:return df
	"""
	for index, row in df.iterrows():
		vec = np.array([df.loc[index,"x_coord"], df.loc[index,"y_coord"], df.loc[index,"z_coord"]])
		cross_prod = np.dot(vec, matrix)
		df.at[index,"x_coord"] = cross_prod[0]
		df.at[index,"y_coord"] = cross_prod[1]
		df.at[index,"z_coord"] = cross_prod[2]
	return 	df

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
