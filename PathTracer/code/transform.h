#pragma once

#ifndef TRANSFORMH
#define TRANSFORMH

#include "matrix4.h"
#include "quaternion.h"

class transform
{
private:
	mutable bool m_dirty;       // Should update or not.
	vec3 m_scale;           // Object's scale.
	quaternion m_rotation;      // Object's rotation.
	vec3 m_translation;     // Object's translation.
	matrix4 m_world;          // Object's model matrix.
	matrix4 m_invWorld;       // Object's normal matrix;

public:
	// Object's local axis.
	static const vec3 LocalForward;
	static const vec3 LocalUp;
	static const vec3 LocalRight;

	// ctor/dtor
	transform();
	~transform() = default;

	// Getter.
	matrix4 toMatrix();
	matrix4 toInvMatrix();

	// Transformation.
	void scale(const vec3 &ds);
	void translate(const vec3 &dt);
	void rotate(const vec3 &axis, float angle);
	void setScale(const vec3 &s);
	void setRotation(const quaternion &r);
	void setTranslation(const vec3 &t);

	// Query object's axis.
	vec3 forward() const;
	vec3 up() const;
	vec3 right() const;

	// Transformation getter.
	vec3 translation() const { return m_translation; }
	quaternion rotation() const { return m_rotation; }
	vec3 scale() const { return m_scale; }
	bool getDirtry() const { return m_dirty; }
};

const vec3 transform::LocalUp(0.0f, 1.0f, 0.0f);
const vec3 transform::LocalForward(0.0f, 0.0f, 1.0f);
const vec3 transform::LocalRight(1.0f, 0.0f, 0.0f);

transform::transform()
	:m_dirty(true), m_scale(1.0, 1.0, 1.0)
{
	m_world.loadIdentity();
	m_invWorld.loadIdentity();
}

matrix4 transform::toMatrix()
{
	if (m_dirty)
	{
		m_dirty = false;
		m_world = m_rotation.toMatrix();
		matrix4 trans, scals;
		scals.setScale(m_scale);
		trans.setTranslation(m_translation);
		m_world = trans * m_world * scals;
		m_invWorld = m_world.getInverseTranspose();
	}
	return m_world;
}

matrix4 transform::toInvMatrix()
{
	if (m_dirty)
	{
		m_dirty = false;
		m_world = m_rotation.toMatrix();
		matrix4 trans, scals;
		scals.setScale(m_scale);
		trans.setTranslation(m_translation);
		m_world = trans * m_world * scals;
		m_invWorld = m_world.getInverseTranspose();
	}
	return m_invWorld;
}

void transform::scale(const vec3 &ds)
{
	m_dirty = true;
	m_scale.x *= ds.x;
	m_scale.y *= ds.y;
	m_scale.z *= ds.z;
}

void transform::translate(const vec3 &dt)
{
	m_dirty = true;
	m_translation += dt;
}

void transform::rotate(const vec3 &axis, float angle)
{
	m_dirty = true;
	quaternion newRot;
	newRot.setRotationAxis(axis, angle);
	m_rotation = newRot * m_rotation;
}

void transform::setScale(const vec3 &s)
{
	m_dirty = true;
	m_scale = s;
}

void transform::setRotation(const quaternion &r)
{
	m_dirty = true;
	m_rotation = r;
}

void transform::setTranslation(const vec3 &t)
{
	m_dirty = true;
	m_translation = t;
}

vec3 transform::forward() const
{
	return m_rotation * LocalForward;
}

vec3 transform::up() const
{
	return m_rotation * LocalUp;
}

vec3 transform::right() const
{
	return m_rotation * LocalRight;
}

#endif