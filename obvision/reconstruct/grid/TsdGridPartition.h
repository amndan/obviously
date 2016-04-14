#ifndef TSDGRIDPARTITION_H
#define TSDGRIDPARTITION_H

#include "obcore/math/linalg/linalg.h"
#include "obvision/reconstruct/grid/SensorPolar2D.h"
#include "obvision/reconstruct/grid/TsdGridComponent.h"

namespace obvious
{

/**
 * Cell structure for TsdGrid
 * @param tsd value of the truncated signed distance function
 * @param weight used to calculate contribution of pushed samples
 */
struct TsdCell
{
  obfloat tsd;
  obfloat weight;
};

/**
 * @class TsdGridPartition
 * @brief Acceleration structure for TsdGrid approach
 * @author Stefan May
 */
class TsdGridPartition : public TsdGridComponent
{
  friend class TsdGrid;
public:

  /**
   * Standard constructor
   * Allocates and initializes space and matrices
   * @param[in] x start index in x-dimension
   * @param[in] y start index in y-dimension
   * @param[in] dimY Number of cells in y-dimension
   * @param[in] dimX Number of cells in x-dimension
   * @param[in] dimY Number of cells in y-dimension
   * @param[in] cellSize Size of cell in meters
   */
  TsdGridPartition(const unsigned int x, const unsigned int y, const unsigned int dimX, const unsigned int dimY, const obfloat cellSize);

  /**
   * Destructor
   */
  ~TsdGridPartition();

  /**
   * TSD access operator
   * @param y y-index
   * @param x x-index
   * @return TSD
   */
  obfloat& operator () (unsigned int y, unsigned int x) const { return _grid[y][x].tsd; }

  /**
   * Initialize partitions
   */
  void init(obfloat maxTruncation);

  /**
   * Initialization indication
   * @return initialization flag
   */
  bool isInitialized() const { return _initialized; }

  /**
   * Emptiness indication
   * @return emptiness flag
   */
  bool isEmpty() const { return (!_initialized && _initWeight > 0.0); }

  /**
   * Get x-index
   * @return x-index
   */
  unsigned int getX() const { return _x; }

  /**
   * Get y-index
   * @return y-index
   */
  unsigned int getY() const { return _y; }

  /**
   * Get homogeneous cell coordinates
   * @return cell coordinate matrix of partition
   */
  Matrix* getCellCoordsHom() const { return _cellCoordsHom; }

  /**
   * Get partition coordinates
   * @return partition coordinates matrix
   */
  Matrix* getPartitionCoords() const { return _partCoords; }

  /**
   * Get width, i.e., number of cells in x-dimension
   * @return width
   */
  unsigned int getWidth() const { return _cellsX; }

  /**
   * Get height, i.e., number of cells in y-dimension
   * @return height
   */
  unsigned int getHeight() const { return _cellsY; }

  /**
   * Get number of cells in partition, i.e., width x height
   * @return number of cells
   */
  unsigned int getSize() const { return _cellsX*_cellsY; }

  /**
   * Add TSD value at certain cell
   * @param x cell x-index
   * @param y cell y-index
   * @param sdf SDF
   * @param weight measurement weight
   */
  void addTsd(const unsigned int x, const unsigned int y, const obfloat sdf, const obfloat weight);
  void addTsdForce(const unsigned int x, const unsigned int y, const obfloat sdf, const obfloat weight);
  void deleteTsdForce(const unsigned int x, const unsigned int y, const obfloat sdf, const obfloat weight);
  obfloat sdfToTsd(obfloat sd);

  /**
   * Increase emptiness of whole partition, i.e., every measurement ray passes through partition
   */
  virtual void increaseEmptiness();

  /**
   * Interpolate bilinear within cell
   * @param x x-index
   * @param y y-index
   * @param dx x-coordinate within cell
   * @param dy y-coordinate within cell
   * @return interpolated TSD
   */
  obfloat interpolateBilinear(int x, int y, obfloat dx, obfloat dy);

private:

  static Matrix* _partCoords;

  TsdCell** _grid;

  obfloat _cellSize;

  Matrix* _cellCoordsHom;

  unsigned int _cellsX;

  unsigned int _cellsY;

  unsigned int _x;

  unsigned int _y;

  obfloat _initWeight;

  bool _initialized;

  obfloat _maxTruncation;

  obfloat _invMaxTruncation;

  obfloat _eps; // todo: comment

};

inline obfloat TsdGridPartition::sdfToTsd(obfloat sd)
{
  obfloat tsd;

  sd = std::abs(sd);

  tsd = min(sd * _invMaxTruncation, TSDINC);

  if (sd < 0.01)
  {
    tsd = -1;
  }

  return tsd;
}

inline void TsdGridPartition::addTsd(const unsigned int x, const unsigned int y, const obfloat sd, const obfloat weight)
{
  // Factor avoids thin objects to be removed when seen from two sides
  // Todo: Find better solution
  // old tsd function
//  if(sd >= -_maxTruncation)
//  {
//    TsdCell* cell = &_grid[y][x];
//
//    obfloat tsd = min(sd * _invMaxTruncation, TSDINC);
//
//    /**
//     *  The following weighting were proposed by:
//     *  E. Bylow, J. Sturm, C. Kerl, F. Kahl, and D. Cremers.
//     *  Real-time camera tracking and 3d reconstruction using signed distance functions.
//     *  In Robotics: Science and Systems Conference (RSS), June 2013.
//     *  SM: ... we did not achieve a noticeable effect with this.
//     */
//    // obfloat span = -_maxTruncation - _eps;
//    // obfloat sigma = 3.0/(span*span);
//    // obfloat w = 1.0;
//    // if(sd <= _eps) w = exp(-_sigma*(sd-_eps)*(sd-_eps));
//
//    // Experimental: Increase weight for cells close to the surface
//    // If we see a thin surface from both sides, this might prevent temporal removement
//    obfloat w = 0.01;  // todo: magic number
//    if(fabs(sd)<_eps) w = 1.0;
//    w *= weight;
//
//    if(isnan(cell->tsd))  // if cell is empty
//    {
//      cell->tsd = tsd;
//      cell->weight += w;  // weight = 1.0
//    }
//    else
//    {
//      //cell->weight = min(cell->weight+TSDINC, TSDGRIDMAXWEIGHT);
//      //cell->tsd   = (cell->tsd * (cell->weight - TSDINC) + tsd) / cell->weight;
//
//      cell->tsd   = (cell->tsd * cell->weight + tsd * w) / (cell->weight + w);  // weight tsd
//      cell->weight = min(cell->weight+w, TSDGRIDMAXWEIGHT);
//    }
//  }

  // new tsd function 03/04/16
  if(sd >= -_maxTruncation)
    {
      TsdCell* cell = &_grid[y][x]; // cell is old value

      obfloat tsd = min(sd * _invMaxTruncation, TSDINC); // tsd is new value --> here increment is calculated

      //new tsd function 03/04/16
      //obfloat tsd = sdfToTsd(sd);

      /**
       *  The following weighting were proposed by:
       *  E. Bylow, J. Sturm, C. Kerl, F. Kahl, and D. Cremers.
       *  Real-time camera tracking and 3d reconstruction using signed distance functions.
       *  In Robotics: Science and Systems Conference (RSS), June 2013.
       *  SM: ... we did not achieve a noticeable effect with this.
       */
      // obfloat span = -_maxTruncation - _eps;
      // obfloat sigma = 3.0/(span*span);
      // obfloat w = 1.0;
      // if(sd <= _eps) w = exp(-_sigma*(sd-_eps)*(sd-_eps));

      // Experimental: Increase weight for cells close to the surface
      // If we see a thin surface from both sides, this might prevent temporal removement
      obfloat w = 0.01;  // todo: magic number
      if(fabs(sd)<_eps) w = 1.0;
      w *= weight;

      if(isnan(cell->tsd))  // if cell is empty
      {
        cell->tsd = tsd;
        cell->weight += w;  // weight = 1.0
      }
      else
      {
        //cell->weight = min(cell->weight+TSDINC, TSDGRIDMAXWEIGHT);
        //cell->tsd   = (cell->tsd * (cell->weight - TSDINC) + tsd) / cell->weight;

        cell->tsd   = (cell->tsd * cell->weight + tsd * w) / (cell->weight + w);  // weight tsd
        cell->weight = min(cell->weight+w, TSDGRIDMAXWEIGHT);
      }
    }

}

inline void TsdGridPartition::addTsdForce(const unsigned int x, const unsigned int y, const obfloat sd, const obfloat weight)
{
  // Factor avoids thin objects to be removed when seen from two sides
  // Todo: Find better solution
  // include the TSD-Function around the forced point
  // old tsd function
//  if(sd >= -_maxTruncation)
//  {
//    TsdCell* cell = &_grid[y][x];
//
//    obfloat tsd = min(sd * _invMaxTruncation, TSDINC);
//
//    cell->tsd = tsd;
//    cell->weight = TSDGRIDMAXWEIGHT;
//  }
//  // defines a unknown area behind the forced point => no second wall.
//  if((sd >= -(1.5*_maxTruncation)) and (sd <= -_maxTruncation))
//  {
//    TsdCell* cell = &_grid[y][x];
//
//    cell->tsd = NAN;
//    cell->weight = _initWeight;
//  }

  // new tsd function 03/04/16
//    if(sd >= -_maxTruncation)
//    {
//      TsdCell* cell = &_grid[y][x];
//
//      obfloat tsd = sdfToTsd(sd);
//
//      cell->tsd = tsd;
//      cell->weight = TSDGRIDMAXWEIGHT;
//    }
//    // defines a unknown area behind the forced point => no second wall.
//    if((sd >= -(1.5*_maxTruncation)) and (sd <= -_maxTruncation))
//    {
//      TsdCell* cell = &_grid[y][x];
//
//      cell->tsd = 1;
//      cell->weight = _initWeight;
//    }

    // new tsd function 04/11/16
      if(sd >= -_maxTruncation)
      {
        TsdCell* cell = &_grid[y][x];

          obfloat tsd = min(sd * _invMaxTruncation, TSDINC); // tsd is new value --> here increment is calculated

        cell->tsd = tsd;
        cell->weight = TSDGRIDMAXWEIGHT;
      }
      // defines a unknown area behind the forced point => no second wall.
      if((sd >= -(2*_maxTruncation)) and (sd < -_maxTruncation))
      {
        TsdCell* cell = &_grid[y][x];
          obfloat tsd = min(sd * _invMaxTruncation, TSDINC); // tsd is new value --> here increment is calculated

        cell->tsd = 1-tsd;
        cell->weight = _initWeight;
      }
}

inline void TsdGridPartition::deleteTsdForce(const unsigned int x, const unsigned int y, const obfloat sd, const obfloat weight)
{
  // Factor avoids thin objects to be removed when seen from two sides
  // Todo: Find better solution
  // erase the point and sets the area of the TSD-Function as unknown
  if((sd >= -_maxTruncation) and (sd <= _maxTruncation))
  {
    TsdCell* cell = &_grid[y][x];

    cell->tsd = NAN;
    cell->weight = _initWeight;
  }


}

inline obfloat TsdGridPartition::interpolateBilinear(int x, int y, obfloat dx, obfloat dy)
{
  // Interpolate
  return _grid[y][x].tsd         * (1. - dy) * (1. - dx)
       + _grid[y + 1][x + 0].tsd *       dy  * (1. - dx)
       + _grid[y + 0][x + 1].tsd * (1. - dy) *       dx
       + _grid[y + 1][x + 1].tsd *       dy  *       dx;
}

}

#endif
