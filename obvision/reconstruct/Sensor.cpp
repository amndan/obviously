#ifndef SENSOR_H
#define SENSOR_H

#include "obcore/math/Matrix.h"
#include "Sensor.h"
#include <string.h>

namespace obvious
{

Sensor::Sensor(unsigned int dim)
{
  _dim = dim;

  _Pose = new Matrix(_dim+1, _dim+1);
  _Pose->setIdentity();
}

Sensor::~Sensor()
{
  delete _Pose;
}

unsigned int Sensor::getRealMeasurementSize()
{
  return _size;
}

void Sensor::setRealMeasurementData(double* data)
{
  memcpy(_data, data, _size*sizeof(*data));
}

double* Sensor::getRealMeasurementData()
{
  return _data;
}

void Sensor::setRealMeasurementMask(bool* mask)
{
  memcpy(_mask, mask, _size*sizeof(*mask));
}

bool* Sensor::getRealMeasurementMask()
{
  return _mask;
}

void Sensor::transform(Matrix* T)
{
  (*_Pose) *= (*T);
}

Matrix* Sensor::getPose()
{
  return _Pose;
}

void Sensor::getPosition(double* tr)
{
  tr[0] = (*_Pose)[0][_dim];
  tr[1] = (*_Pose)[1][_dim];
  if(_dim==3)
    tr[2] = (*_Pose)[2][_dim];
}

}

#endif
