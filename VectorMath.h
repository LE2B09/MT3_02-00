#pragma once

#include "Matrix4x4.h"
#include "Vector3.h"
#include <assert.h>
#include <cmath>

// 加算
static Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result{};
	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;
	return result;
}

// 減算
static Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result{};
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;
	return result;
}

// スカラー倍
static Vector3 Multiply(float scalar, const Vector3& v) {
	Vector3 result{};
	result.x = scalar * v.x;
	result.y = scalar * v.y;
	result.z = scalar * v.z;
	return result;
}

// 内積
static float Dot(const Vector3& v1, const Vector3& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }

// 長さ（ノルム）
static float Length(const Vector3& v) { return sqrtf(powf(v.x, 2) + powf(v.y, 2) + powf(v.z, 2)); }

//// 正規化
//static Vector3 Nomalize(const Vector3& v) {
//	float length = Length(v);
//	Vector3 result{};
//	if (length != 0.0) {
//		result.x = v.x / length;
//		result.y = v.y / length;
//		result.z = v.z / length;
//	}
//	return result;
//}

// 座標変換
static Vector3 Transforms(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result{};
	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];
	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;
	return result;
}

//// クロス積
//static Vector3 Cross(const Vector3& v1, const Vector3& v2) {
//	Vector3 result{};
//	result.x = v1.y * v2.z - v1.z * v2.y;
//	result.y = v1.z * v2.x - v1.x * v2.z;
//	result.z = v1.x * v2.y - v1.y * v2.x;
//	return result;
//}

//正射影ベクトル（ベクトル射影）
static Vector3 Project(const Vector3& v1, const Vector3& v2)
{
	float lengthSqr = Length(v2);
	if (lengthSqr == 0.0f)
	{
		return Vector3{ 0.0f,0.0f,0.0f };
	}

	return Multiply(Dot(v1, v2) / powf(lengthSqr, 2), v2);
}

////最近接点
//static Vector3 ClosestPoint(const Vector3& point, const Segment& segment)
//{
//	return Subtract(segment.origin, Project(point, segment.diff));
//}
//
////最近接点を計算し、クランプ処理を行う関数
//static Vector3 ClosestPointWithClamp(const Vector3& point, const Segment& segment)
//{
//	//最近接点を計算
//	Vector3 closestPoint = ClosestPoint(point, segment);
//
//	//線分の始点から終点へのベクトル
//	Vector3 segmentVector = segment.diff;
//
//	//線分の始点から最近接点へのベクトル
//	Vector3 closestPointVector = Subtract(closestPoint, segment.origin);
//
//	//線分の始点から終点への長さの二乗
//	float segmentLengthSquared = Length(segmentVector);
//
//	//ベクトルの内積を計算
//	float dotProduct = Dot(segmentVector, closestPointVector);
//
//	//最近接点が線分の始点よりも前にある場合、線分の始点を最近接点に設定
//	if (dotProduct < 0.0f)
//	{
//		closestPoint = segment.origin;
//	}
//
//	//最近接点が線分の始点よりも後ろにある場合、線分の終点を最近接点に設定
//	else if (dotProduct > segmentLengthSquared)
//	{
//		closestPoint = Add(segment.origin, segmentVector);
//	}
//
//	return closestPoint;
//}
