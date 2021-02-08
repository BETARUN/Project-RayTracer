#pragma once

#ifndef QUATERNIONH
#define QUATERNIONH

#include "vec3.h"
#include "matrix4.h"

class quaternion
{

public:
	const static quaternion identity;
	double x, y, z, w;

	quaternion();
	quaternion(float x, float y, float z, float w);
	quaternion(float yaw, float pitch, float roll);
	~quaternion() = default;

	void set(float _x, float _y, float _z, float _w);
	void setEulerAngle(float yaw, float pitch, float roll);
	void setRotationAxis(vec3 axis, double angle);

	quaternion inverse() const;
	quaternion conjugate() const;
	vec3 eulerAngle() const;
	matrix4 toMatrix() const;

	static float dot(const quaternion &lhs, const quaternion &rhs);
	static quaternion lerp(const quaternion &a, const quaternion &b, float t);
	static quaternion slerp(const quaternion &a, const quaternion &b, float t);
	static float angle(const quaternion &lhs, const quaternion &rhs);

	void operator*(float s);
	void operator+(const quaternion &q);
	void operator-(const quaternion &q);

	friend quaternion operator * (const quaternion& lhs, const quaternion& rhs);
	friend vec3 operator *(const quaternion& rotation, const vec3& point);

private:
	vec3 eulerAngles;
};

const quaternion quaternion::identity(0, 0, 0, 1);

quaternion::quaternion()
{
	x = y = z = 0;
	w = 1;
}

quaternion::quaternion(float _x, float _y, float _z, float _w)
{
	double mag = _x * _x + _y * _y + _z * _z + _w * _w;
	x = _x / mag;
	y = _y / mag;
	z = _z / mag;
	w = _w / mag;
}

quaternion::quaternion(float yaw, float pitch, float roll)
{
	this->setEulerAngle(yaw, pitch, roll);
}

void quaternion::setEulerAngle(float yaw, float pitch, float roll)
{
	float  angle;
	float  sinRoll, sinPitch, sinYaw, cosRoll, cosPitch, cosYaw;

	angle = yaw * 0.5f;
	sinYaw = sin(angle);
	cosYaw = cos(angle);

	angle = pitch * 0.5f;
	sinPitch = sin(angle);
	cosPitch = cos(angle);

	angle = roll * 0.5f;
	sinRoll = sin(angle);
	cosRoll = cos(angle);

	float _x = cosRoll * cosPitch*sinYaw - sinRoll * sinPitch*cosYaw;
	float _y = cosRoll * sinPitch*cosYaw + sinRoll * cosPitch*sinYaw;
	float _z = sinRoll * cosPitch*cosYaw - cosRoll * sinPitch*sinYaw;
	float _w = cosRoll * cosPitch*cosYaw + sinRoll * sinPitch*sinYaw;

	float mag = _x * _x + _y * _y + _z * _z + _w * _w;
	x = _x / mag;
	y = _y / mag;
	z = _z / mag;
	w = _w / mag;
}

void quaternion::setRotationAxis(vec3 axis, double angle)
{
	angle = radians(angle);
	axis.make_unit_vector();
	double angleDiv2 = angle * 0.5;
	double sinAngle = sin(angleDiv2);
	x = axis.x * sinAngle;
	y = axis.y * sinAngle;
	z = axis.z * sinAngle;
	w = cos(angleDiv2);
}

quaternion quaternion::inverse() const
{
	return quaternion(-x, -y, -z, w);
}

quaternion quaternion::conjugate() const
{
	return quaternion(-x, -y, -z, w);
}

vec3 quaternion::eulerAngle() const
{
	float yaw = atan2(2 * (w * x + z * y), 1 - 2 * (x * x + y * y));
	float pitch = asin(2 * (w * y - x * z));
	float roll = atan2(2 * (w * z + x * y), 1 - 2 * (z * z + y * y));
	//    if(pitch < -1.0f)pitch = -1.0f;
	//    if(pitch > +1.0f)pitch = +1.0f;
	return vec3(
		angles(yaw),
		angles(pitch),
		angles(roll));
}

matrix4 quaternion::toMatrix() const
{
	matrix4 result(
		1.0f - 2.0f*y*y - 2.0f*z*z, 2.0f*x*y - 2.0f*z*w, 2.0f*x*z + 2.0f*y*w, 0.0f,
		2.0f*x*y + 2.0f*z*w, 1.0f - 2.0f*x*x - 2.0f*z*z, 2.0f*y*z - 2.0f*x*w, 0.0f,
		2.0f*x*z - 2.0f*y*w, 2.0f*y*z + 2.0f*x*w, 1.0f - 2.0f*x*x - 2.0f*y*y, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
	result.transpose();
	return result;
}

void quaternion::set(float _x, float _y, float _z, float _w)
{
	x = _x; y = _y; z = _z; w = _w;
}

float quaternion::dot(const quaternion &lhs, const quaternion &rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
}

quaternion quaternion::lerp(const quaternion &a, const quaternion &b, float t)
{
	return quaternion(
		(1 - t) * a.x + t * b.x,
		(1 - t) * a.y + t * b.y,
		(1 - t) * a.z + t * b.z,
		(1 - t) * a.w + t * b.w
	);
}

quaternion quaternion::slerp(const quaternion &a, const quaternion &b, float t)
{
	float cos_theta = dot(a, b);

	// if B is on opposite hemisphere from A, use -B instead
	float sign;
	if (cos_theta < 0.f)
	{
		cos_theta = -cos_theta;
		sign = -1.f;
	}
	else sign = 1.f;

	float c1, c2;
	if (cos_theta > 1.f - 0.000001f)
	{
		// if q2 is (within precision limits) the same as q1,
		// just linear interpolate between A and B
		c2 = t;
		c1 = 1.f - t;
	}
	else
	{
		float theta = acos(cos_theta);
		float sin_theta = sin(theta);
		float t_theta = t * theta;
		float inv_sin_theta = 1.f / sin_theta;
		c2 = sin(t_theta) * inv_sin_theta;
		c1 = sin(theta - t_theta) * inv_sin_theta;
	}
	c2 *= sign;
	// interpolate
	return quaternion(
		a.x * c1 + b.x * c2,
		a.y * c1 + b.y * c2,
		a.z * c1 + b.z * c2,
		a.w * c1 + b.w * c2);
}

float quaternion::angle(const quaternion &lhs, const quaternion &rhs)
{
	float cos_theta = dot(lhs, rhs);
	if (cos_theta < 0.0f)
		cos_theta = -cos_theta;
	return 2.0f * angles(acos(cos_theta));
}

void quaternion::operator*(float s)
{
	x *= s; y *= s; z *= s; w *= s;
}

void quaternion::operator+(const quaternion &q)
{
	x += q.x; y += q.y; z += q.z; w += q.w;
}

void quaternion::operator-(const quaternion &q)
{
	x -= q.x; y -= q.y; z -= q.z; w -= q.w;
}

quaternion operator *(const quaternion &lhs, const quaternion &rhs)
{
	float w1 = lhs.w;
	float w2 = rhs.w;
	vec3 v1(lhs.x, lhs.y, lhs.z);
	vec3 v2(rhs.x, rhs.y, rhs.z);
	float w3 = w1 * w2 - dot(v1, v2);
	vec3 v3 = cross(v1, v2) + v2 * w1 + v1 * w2;
	return quaternion(v3.x, v3.y, v3.z, w3);
}

vec3 operator *(const quaternion &q, const vec3 &v)
{
	// Extract the vector part of the quaternion
	vec3 u(q.x, q.y, q.z);
	// Extract the scalar part of the quaternion
	float s = q.w;
	// Do the math
	return u * 2.0f * dot(u, v)
		+ v * (s*s - dot(u, u))
		+ cross(u, v) * 2.0f * s;
}

#endif
