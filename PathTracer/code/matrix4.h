#pragma once
#ifndef matrix4H
#define matrix4H

#include "mathutils.h"
#include "vec4.h"

class matrix4
{
public:
	float entries[16];

	// constructors
	matrix4() { loadIdentity(); }
	matrix4(float e0, float e1, float e2, float e3,
		float e4, float e5, float e6, float e7,
		float e8, float e9, float e10, float e11,
		float e12, float e13, float e14, float e15);
	matrix4(const float *rhs);
	matrix4(const matrix4 &rhs);
	~matrix4() = default;

	// setter,getter
	void setEntry(int position, float value);
	float getEntry(int position) const;
	vec4 getRow(int position) const;
	vec4 getColumn(int position) const;
	void loadIdentity();
	void loadZero();

	// overloaded operators
	matrix4 operator+(const matrix4 & rhs) const;
	matrix4 operator-(const matrix4 & rhs) const;
	matrix4 operator*(const matrix4 & rhs) const;
	matrix4 operator*(const float rhs) const;
	matrix4 operator/(const float rhs) const;

	bool operator==(const matrix4 & rhs) const;
	bool operator!=(const matrix4 & rhs) const;

	void operator+=(const matrix4 & rhs);
	void operator-=(const matrix4 & rhs);
	void operator*=(const matrix4 & rhs);
	void operator*=(const float rhs);
	void operator/=(const float rhs);

	matrix4 operator-() const;
	matrix4 operator+() const { return (*this); }
	vec4 operator*(const vec4 rhs) const;

	// inverse, transpose
	void inverted();
	matrix4 getInverse() const;
	void transpose();
	matrix4 getTranspose() const;
	void invertTranspose();
	matrix4 getInverseTranspose() const;

	// operation on space
	void setTranslation(const vec3 & translation);
	void setScale(const vec3 & scaleFactor);
	void setRotationAxis(const double angle, const vec3 & axis);
	void setRotationX(const double angle);
	void setRotationY(const double angle);
	void setRotationZ(const double angle);
	void setRotationEuler(const double angleX, const double angleY, const double angleZ);
	void setPerspective(float fovy, float aspect, float near, float far);
	void setOrtho(float left, float right, float bottom, float top, float near, float far);
	void setLookAt(vec3 cameraPos, vec3 target, vec3 worldUp);
	void setViewPort(int left, int bottom, int width, int height);
};

matrix4::matrix4(float e0, float e1, float e2, float e3,
	float e4, float e5, float e6, float e7,
	float e8, float e9, float e10, float e11,
	float e12, float e13, float e14, float e15)
{
	entries[0] = e0;
	entries[1] = e1;
	entries[2] = e2;
	entries[3] = e3;
	entries[4] = e4;
	entries[5] = e5;
	entries[6] = e6;
	entries[7] = e7;
	entries[8] = e8;
	entries[9] = e9;
	entries[10] = e10;
	entries[11] = e11;
	entries[12] = e12;
	entries[13] = e13;
	entries[14] = e14;
	entries[15] = e15;
}

matrix4::matrix4(const float *rhs)
{
	memcpy(entries, rhs, 16 * sizeof(float));
}

matrix4::matrix4(const matrix4 &rhs)
{
	memcpy(entries, rhs.entries, 16 * sizeof(float));
}

void matrix4::setEntry(int position, float value)
{
	if (position >= 0 && position <= 15)
		entries[position] = value;
}

float matrix4::getEntry(int position) const
{
	if (position >= 0 && position <= 15)
		return entries[position];
	return 0.0f;
}

vec4 matrix4::getRow(int position) const
{
	switch (position)
	{
	case 0:
		return vec4(entries[0], entries[4], entries[8], entries[12]);
	case 1:
		return vec4(entries[1], entries[5], entries[9], entries[13]);
	case 2:
		return vec4(entries[2], entries[6], entries[10], entries[14]);
	case 3:
		return vec4(entries[3], entries[7], entries[11], entries[15]);
	}
	return vec4(0.0f, 0.0f, 0.0f, 0.0f);
}

vec4 matrix4::getColumn(int position) const
{
	switch (position)
	{
	case 0:
		return vec4(entries[0], entries[1], entries[2], entries[3]);
	case 1:
		return vec4(entries[4], entries[5], entries[6], entries[7]);
	case 2:
		return vec4(entries[8], entries[9], entries[10], entries[11]);
	case 3:
		return vec4(entries[12], entries[13], entries[14], entries[15]);
	}
	return vec4(0.0f, 0.0f, 0.0f, 0.0f);
}

void matrix4::loadIdentity()
{
	memset(entries, 0, 16 * sizeof(float));
	entries[0] = 1.0f;
	entries[5] = 1.0f;
	entries[10] = 1.0f;
	entries[15] = 1.0f;
}

void matrix4::loadZero()
{
	memset(entries, 0, 16 * sizeof(float));
}

matrix4 matrix4::operator+(const matrix4 &rhs) const
{
	return matrix4
	(entries[0] + rhs.entries[0],
		entries[1] + rhs.entries[1],
		entries[2] + rhs.entries[2],
		entries[3] + rhs.entries[3],
		entries[4] + rhs.entries[4],
		entries[5] + rhs.entries[5],
		entries[6] + rhs.entries[6],
		entries[7] + rhs.entries[7],
		entries[8] + rhs.entries[8],
		entries[9] + rhs.entries[9],
		entries[10] + rhs.entries[10],
		entries[11] + rhs.entries[11],
		entries[12] + rhs.entries[12],
		entries[13] + rhs.entries[13],
		entries[14] + rhs.entries[14],
		entries[15] + rhs.entries[15]);
}

matrix4 matrix4::operator-(const matrix4 &rhs) const
{
	return matrix4
	(entries[0] - rhs.entries[0],
		entries[1] - rhs.entries[1],
		entries[2] - rhs.entries[2],
		entries[3] - rhs.entries[3],
		entries[4] - rhs.entries[4],
		entries[5] - rhs.entries[5],
		entries[6] - rhs.entries[6],
		entries[7] - rhs.entries[7],
		entries[8] - rhs.entries[8],
		entries[9] - rhs.entries[9],
		entries[10] - rhs.entries[10],
		entries[11] - rhs.entries[11],
		entries[12] - rhs.entries[12],
		entries[13] - rhs.entries[13],
		entries[14] - rhs.entries[14],
		entries[15] - rhs.entries[15]);
}

matrix4 matrix4::operator*(const matrix4 &rhs) const
{
	return matrix4
	(entries[0] * rhs.entries[0] + entries[4] * rhs.entries[1] + entries[8] * rhs.entries[2] + entries[12] * rhs.entries[3],
		entries[1] * rhs.entries[0] + entries[5] * rhs.entries[1] + entries[9] * rhs.entries[2] + entries[13] * rhs.entries[3],
		entries[2] * rhs.entries[0] + entries[6] * rhs.entries[1] + entries[10] * rhs.entries[2] + entries[14] * rhs.entries[3],
		entries[3] * rhs.entries[0] + entries[7] * rhs.entries[1] + entries[11] * rhs.entries[2] + entries[15] * rhs.entries[3],
		entries[0] * rhs.entries[4] + entries[4] * rhs.entries[5] + entries[8] * rhs.entries[6] + entries[12] * rhs.entries[7],
		entries[1] * rhs.entries[4] + entries[5] * rhs.entries[5] + entries[9] * rhs.entries[6] + entries[13] * rhs.entries[7],
		entries[2] * rhs.entries[4] + entries[6] * rhs.entries[5] + entries[10] * rhs.entries[6] + entries[14] * rhs.entries[7],
		entries[3] * rhs.entries[4] + entries[7] * rhs.entries[5] + entries[11] * rhs.entries[6] + entries[15] * rhs.entries[7],
		entries[0] * rhs.entries[8] + entries[4] * rhs.entries[9] + entries[8] * rhs.entries[10] + entries[12] * rhs.entries[11],
		entries[1] * rhs.entries[8] + entries[5] * rhs.entries[9] + entries[9] * rhs.entries[10] + entries[13] * rhs.entries[11],
		entries[2] * rhs.entries[8] + entries[6] * rhs.entries[9] + entries[10] * rhs.entries[10] + entries[14] * rhs.entries[11],
		entries[3] * rhs.entries[8] + entries[7] * rhs.entries[9] + entries[11] * rhs.entries[10] + entries[15] * rhs.entries[11],
		entries[0] * rhs.entries[12] + entries[4] * rhs.entries[13] + entries[8] * rhs.entries[14] + entries[12] * rhs.entries[15],
		entries[1] * rhs.entries[12] + entries[5] * rhs.entries[13] + entries[9] * rhs.entries[14] + entries[13] * rhs.entries[15],
		entries[2] * rhs.entries[12] + entries[6] * rhs.entries[13] + entries[10] * rhs.entries[14] + entries[14] * rhs.entries[15],
		entries[3] * rhs.entries[12] + entries[7] * rhs.entries[13] + entries[11] * rhs.entries[14] + entries[15] * rhs.entries[15]);
}

matrix4 matrix4::operator*(const float rhs) const
{
	return matrix4
	(entries[0] * rhs,
		entries[1] * rhs,
		entries[2] * rhs,
		entries[3] * rhs,
		entries[4] * rhs,
		entries[5] * rhs,
		entries[6] * rhs,
		entries[7] * rhs,
		entries[8] * rhs,
		entries[9] * rhs,
		entries[10] * rhs,
		entries[11] * rhs,
		entries[12] * rhs,
		entries[13] * rhs,
		entries[14] * rhs,
		entries[15] * rhs);
}

matrix4 matrix4::operator/(const float rhs) const
{
	if (equal(rhs, 0.0f) || equal(rhs, 1.0f))
		return (*this);
	float temp = 1 / rhs;
	return (*this)*temp;
}

bool matrix4::operator==(const matrix4 &rhs) const
{
	for (int i = 0; i < 16; i++)
	{
		if (!equal(entries[i], rhs.entries[i]))
			return false;
	}
	return true;
}

bool matrix4::operator!=(const matrix4 &rhs) const
{
	return !((*this) == rhs);
}

void matrix4::operator+=(const matrix4 &rhs)
{
	(*this) = (*this) + rhs;
}

void matrix4::operator-=(const matrix4 &rhs)
{
	(*this) = (*this) - rhs;
}

void matrix4::operator*=(const matrix4 &rhs)
{
	(*this) = (*this)*rhs;
}

void matrix4::operator*=(const float rhs)
{
	(*this) = (*this)*rhs;
}

void matrix4::operator/=(const float rhs)
{
	(*this) = (*this) / rhs;
}

matrix4 matrix4::operator-() const
{
	matrix4 result(*this);
	for (int i = 0; i < 16; i++)
		result.entries[i] = -result.entries[i];
	return result;
}

vec4 matrix4::operator*(const vec4 rhs) const
{
	return vec4(
		entries[0] * rhs.x
		+ entries[4] * rhs.y
		+ entries[8] * rhs.z
		+ entries[12] * rhs.w,

		entries[1] * rhs.x
		+ entries[5] * rhs.y
		+ entries[9] * rhs.z
		+ entries[13] * rhs.w,

		entries[2] * rhs.x
		+ entries[6] * rhs.y
		+ entries[10] * rhs.z
		+ entries[14] * rhs.w,

		entries[3] * rhs.x
		+ entries[7] * rhs.y
		+ entries[11] * rhs.z
		+ entries[15] * rhs.w);
}

void matrix4::inverted()
{
	*this = getInverse();
}

matrix4 matrix4::getInverse() const
{
	matrix4 result = getInverseTranspose();
	result.transpose();
	return result;
}

void matrix4::transpose()
{
	*this = getTranspose();
}

matrix4 matrix4::getTranspose() const
{
	return matrix4
	(entries[0], entries[4], entries[8], entries[12],
		entries[1], entries[5], entries[9], entries[13],
		entries[2], entries[6], entries[10], entries[14],
		entries[3], entries[7], entries[11], entries[15]);
}

void matrix4::invertTranspose()
{
	*this = getInverseTranspose();
}

matrix4 matrix4::getInverseTranspose() const
{
	matrix4 result;
	float tmp[12];
	// determinant
	float det;

	// calculate pairs for first 8 elements (cofactors)
	tmp[0] = entries[10] * entries[15];
	tmp[1] = entries[11] * entries[14];
	tmp[2] = entries[9] * entries[15];
	tmp[3] = entries[11] * entries[13];
	tmp[4] = entries[9] * entries[14];
	tmp[5] = entries[10] * entries[13];
	tmp[6] = entries[8] * entries[15];
	tmp[7] = entries[11] * entries[12];
	tmp[8] = entries[8] * entries[14];
	tmp[9] = entries[10] * entries[12];
	tmp[10] = entries[8] * entries[13];
	tmp[11] = entries[9] * entries[12];

	// calculate first 8 elements (cofactors)
	result.setEntry(0, tmp[0] * entries[5] + tmp[3] * entries[6] + tmp[4] * entries[7]
		- tmp[1] * entries[5] - tmp[2] * entries[6] - tmp[5] * entries[7]);

	result.setEntry(1, tmp[1] * entries[4] + tmp[6] * entries[6] + tmp[9] * entries[7]
		- tmp[0] * entries[4] - tmp[7] * entries[6] - tmp[8] * entries[7]);

	result.setEntry(2, tmp[2] * entries[4] + tmp[7] * entries[5] + tmp[10] * entries[7]
		- tmp[3] * entries[4] - tmp[6] * entries[5] - tmp[11] * entries[7]);

	result.setEntry(3, tmp[5] * entries[4] + tmp[8] * entries[5] + tmp[11] * entries[6]
		- tmp[4] * entries[4] - tmp[9] * entries[5] - tmp[10] * entries[6]);

	result.setEntry(4, tmp[1] * entries[1] + tmp[2] * entries[2] + tmp[5] * entries[3]
		- tmp[0] * entries[1] - tmp[3] * entries[2] - tmp[4] * entries[3]);

	result.setEntry(5, tmp[0] * entries[0] + tmp[7] * entries[2] + tmp[8] * entries[3]
		- tmp[1] * entries[0] - tmp[6] * entries[2] - tmp[9] * entries[3]);

	result.setEntry(6, tmp[3] * entries[0] + tmp[6] * entries[1] + tmp[11] * entries[3]
		- tmp[2] * entries[0] - tmp[7] * entries[1] - tmp[10] * entries[3]);

	result.setEntry(7, tmp[4] * entries[0] + tmp[9] * entries[1] + tmp[10] * entries[2]
		- tmp[5] * entries[0] - tmp[8] * entries[1] - tmp[11] * entries[2]);

	//calculate pairs for second 8 elements (cofactors)
	tmp[0] = entries[2] * entries[7];
	tmp[1] = entries[3] * entries[6];
	tmp[2] = entries[1] * entries[7];
	tmp[3] = entries[3] * entries[5];
	tmp[4] = entries[1] * entries[6];
	tmp[5] = entries[2] * entries[5];
	tmp[6] = entries[0] * entries[7];
	tmp[7] = entries[3] * entries[4];
	tmp[8] = entries[0] * entries[6];
	tmp[9] = entries[2] * entries[4];
	tmp[10] = entries[0] * entries[5];
	tmp[11] = entries[1] * entries[4];

	//calculate second 8 elements (cofactors)
	result.setEntry(8, tmp[0] * entries[13] + tmp[3] * entries[14] + tmp[4] * entries[15]
		- tmp[1] * entries[13] - tmp[2] * entries[14] - tmp[5] * entries[15]);

	result.setEntry(9, tmp[1] * entries[12] + tmp[6] * entries[14] + tmp[9] * entries[15]
		- tmp[0] * entries[12] - tmp[7] * entries[14] - tmp[8] * entries[15]);

	result.setEntry(10, tmp[2] * entries[12] + tmp[7] * entries[13] + tmp[10] * entries[15]
		- tmp[3] * entries[12] - tmp[6] * entries[13] - tmp[11] * entries[15]);

	result.setEntry(11, tmp[5] * entries[12] + tmp[8] * entries[13] + tmp[11] * entries[14]
		- tmp[4] * entries[12] - tmp[9] * entries[13] - tmp[10] * entries[14]);

	result.setEntry(12, tmp[2] * entries[10] + tmp[5] * entries[11] + tmp[1] * entries[9]
		- tmp[4] * entries[11] - tmp[0] * entries[9] - tmp[3] * entries[10]);

	result.setEntry(13, tmp[8] * entries[11] + tmp[0] * entries[8] + tmp[7] * entries[10]
		- tmp[6] * entries[10] - tmp[9] * entries[11] - tmp[1] * entries[8]);

	result.setEntry(14, tmp[6] * entries[9] + tmp[11] * entries[11] + tmp[3] * entries[8]
		- tmp[10] * entries[11] - tmp[2] * entries[8] - tmp[7] * entries[9]);

	result.setEntry(15, tmp[10] * entries[10] + tmp[4] * entries[8] + tmp[9] * entries[9]
		- tmp[8] * entries[9] - tmp[11] * entries[10] - tmp[5] * entries[8]);

	// calculate determinant
	det = entries[0] * result.getEntry(0)
		+ entries[1] * result.getEntry(1)
		+ entries[2] * result.getEntry(2)
		+ entries[3] * result.getEntry(3);

	if (det == 0.0f)
		return matrix4();

	result = result / det;
	return result;
}

void matrix4::setTranslation(const vec3 &translation)
{
	loadIdentity();
	entries[12] = translation.x;
	entries[13] = translation.y;
	entries[14] = translation.z;
}

void matrix4::setScale(const vec3 &scaleFactor)
{
	loadIdentity();
	entries[0] = scaleFactor.x;
	entries[5] = scaleFactor.y;
	entries[10] = scaleFactor.z;
}

void matrix4::setRotationAxis(const double angle, const vec3 &axis)
{
	vec3 u = unit_vector(axis);

	float sinAngle = static_cast<float>(sin(M_PI*angle / 180));
	float cosAngle = static_cast<float>(cos(M_PI*angle / 180));
	float oneMinusCosAngle = 1.0f - cosAngle;

	loadIdentity();

	entries[0] = (u.x)*(u.x) + cosAngle * (1 - (u.x)*(u.x));
	entries[4] = (u.x)*(u.y)*(oneMinusCosAngle)-sinAngle * u.z;
	entries[8] = (u.x)*(u.z)*(oneMinusCosAngle)+sinAngle * u.y;

	entries[1] = (u.x)*(u.y)*(oneMinusCosAngle)+sinAngle * u.z;
	entries[5] = (u.y)*(u.y) + cosAngle * (1 - (u.y)*(u.y));
	entries[9] = (u.y)*(u.z)*(oneMinusCosAngle)-sinAngle * u.x;

	entries[2] = (u.x)*(u.z)*(oneMinusCosAngle)-sinAngle * u.y;
	entries[6] = (u.y)*(u.z)*(oneMinusCosAngle)+sinAngle * u.x;
	entries[10] = (u.z)*(u.z) + cosAngle * (1 - (u.z)*(u.z));
}

void matrix4::setRotationX(const double angle)
{
	loadIdentity();
	entries[5] = static_cast<float>(cos(M_PI*angle / 180));
	entries[6] = static_cast<float>(sin(M_PI*angle / 180));
	entries[9] = -entries[6];
	entries[10] = entries[5];
}

void matrix4::setRotationY(const double angle)
{
	loadIdentity();
	entries[0] = static_cast<float>(cos(M_PI*angle / 180));
	entries[2] = -static_cast<float>(sin(M_PI*angle / 180));
	entries[8] = -entries[2];
	entries[10] = entries[0];
}

void matrix4::setRotationZ(const double angle)
{
	loadIdentity();
	entries[0] = static_cast<float>(cos(M_PI*angle / 180));
	entries[1] = static_cast<float>(sin(M_PI*angle / 180));
	entries[4] = -entries[1];
	entries[5] = entries[0];
}

void matrix4::setRotationEuler(const double angleX, const double angleY, const double angleZ)
{
	// pitch, raw, roll
	loadIdentity();
	double cr = cos(M_PI*angleX / 180);
	double sr = sin(M_PI*angleX / 180);
	double cp = cos(M_PI*angleY / 180);
	double sp = sin(M_PI*angleY / 180);
	double cy = cos(M_PI*angleZ / 180);
	double sy = sin(M_PI*angleZ / 180);

	entries[0] = static_cast<float>(cp*cy);
	entries[1] = static_cast<float>(cp*sy);
	entries[2] = static_cast<float>(-sp);

	double srsp = sr * sp;
	double crsp = cr * sp;

	entries[4] = static_cast<float>(srsp*cy - cr * sy);
	entries[5] = static_cast<float>(srsp*sy + cr * cy);
	entries[6] = static_cast<float>(sr*cp);
	entries[8] = static_cast<float>(crsp*cy + sr * sy);
	entries[9] = static_cast<float>(crsp*sy - sr * cy);
	entries[10] = static_cast<float>(cr*cp);
}

void matrix4::setPerspective(float fovy, float aspect, float near, float far)
{
	loadZero();
	// convert fov from degrees to radians
	float rFovy = fovy * M_PI / 180;
	const float tanHalfFovy = tanf(static_cast<float>(rFovy*0.5f));
	entries[0] = 1.0f / (aspect*tanHalfFovy);
	entries[5] = 1.0f / (tanHalfFovy);
	entries[10] = -(far + near) / (far - near);
	entries[11] = -1.0f;
	entries[14] = (-2.0f*near*far) / (far - near);
}


void matrix4::setOrtho(float left, float right, float bottom, float top, float near, float far)
{
	loadIdentity();
	entries[0] = 2.0f / (right - left);
	entries[5] = 2.0f / (top - bottom);
	entries[10] = -2.0f / (far - near);
	entries[12] = -(right + left) / (right - left);
	entries[13] = -(top + bottom) / (top - bottom);
	entries[14] = -(far + near) / (far - near);
}

void matrix4::setLookAt(vec3 cameraPos, vec3 target, vec3 worldUp)
{
	vec3 zAxis = cameraPos - target;
	zAxis.make_unit_vector();
	vec3 xAxis = cross(worldUp,zAxis);
	xAxis.make_unit_vector();
	vec3 yAxis = cross(zAxis, xAxis);
	yAxis.make_unit_vector();

	loadIdentity();
	entries[0] = xAxis.x;
	entries[4] = xAxis.y;
	entries[8] = xAxis.z;

	entries[1] = yAxis.x;
	entries[5] = yAxis.y;
	entries[9] = yAxis.z;

	entries[2] = zAxis.x;
	entries[6] = zAxis.y;
	entries[10] = zAxis.z;

	entries[12] = -dot(xAxis, cameraPos);
	entries[13] = -dot(yAxis, cameraPos);
	entries[14] = -dot(zAxis, cameraPos);
}

void matrix4::setViewPort(int left, int top, int width, int height)
{
	loadIdentity();
	entries[0] = static_cast<float>(width) / 2.0f;
	entries[5] = -static_cast<float>(height) / 2.0f;
	entries[12] = static_cast<float>(left) + static_cast<float>(width) / 2.0f;
	entries[13] = static_cast<float>(top) + static_cast<float>(height) / 2.0f;
}
#endif // Matrix4x4_H