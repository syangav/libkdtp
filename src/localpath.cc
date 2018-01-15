/*
 * Copyright (c) 2016-2018 LAAS/CNRS
 * All rights reserved.
 *
 * Redistribution and use  in source  and binary  forms,  with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *   1. Redistributions of  source  code must retain the  above copyright
 *      notice and this list of conditions.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice and  this list of  conditions in the  documentation and/or
 *      other materials provided with the distribution.
 *
 *                                     Alexandre Boeuf on Fri Aug 19 2016
 */
#include "kdtp.h"

namespace kdtp {

  LocalPath::LocalPath(const Robot &robot,
                       const State &init, const State &end, double duration):
    robot_(robot)
  {
    unsigned int nbDof = robot_.getNbDof();
    double init_[3], end_[3];

    duration_ = duration;
    for(unsigned int i=0; i<nbDof; i++) {
      init_[0] = init.position()[i];
      init_[1] = init.velocity()[i];
      init_[2] = init.acceleration()[i];

      end_[0] = end.position()[i];
      end_[1] = end.velocity()[i];
      end_[2] = end.acceleration()[i];

      splines_.push_back(Spline(robot_.getDof(i), init_, end_));

      if (duration_ < splines_[i].getDuration())
        duration_ = splines_[i].getDuration();
    }

    for(unsigned int i = 0; i < splines_.size(); i++)
      splines_[i].synchronize(duration_);
  }

  std::vector<double>
  LocalPath::getPositionAt(double time) const
  {
    std::vector<double> position(splines_.size());

    for(unsigned int i=0; i < splines_.size(); i++)
      position[i] = splines_[i].getPositionAt(time);

    return position;
  }

  std::vector<double>
  LocalPath::getVelocityAt(double time) const
  {
    std::vector<double> velocity(splines_.size());

    for(unsigned int i=0; i < splines_.size(); i++)
      velocity[i] = splines_[i].getVelocityAt(time);

    return velocity;
  }

  std::vector<double>
  LocalPath::getAccelerationAt(double time) const
  {
    std::vector<double> acceleration(splines_.size());

    for(unsigned int i=0; i< splines_.size(); i++)
      acceleration[i] = splines_[i].getAccelerationAt(time);

    return acceleration;
  }

  std::vector<double>
  LocalPath::getJerkAt(double time) const
  {
    std::vector<double> jerk(splines_.size());

    for(unsigned int i=0; i < splines_.size(); i++)
      jerk[i] = splines_[i].getJerkAt(time);

    return jerk;
  }

  std::vector<double>
  LocalPath::getSnapAt(double time) const
  {
    std::vector<double> snap(splines_.size());

    for(unsigned int i=0; i < splines_.size(); i++)
      snap[i] = splines_[i].getSnapAt(time);

    return snap;
  }

  std::vector<std::vector<double> >
  LocalPath::getAllAt(double time) const
  {
    std::vector<std::vector<double> > all(splines_.size());

    for(unsigned int i=0; i < splines_.size(); i++)
      all[i] = splines_[i].getAllAt(time);

    return all;
  }

  State
  LocalPath::getStateAt(double time) const
  {
    State state(robot_);

    for(unsigned int i=0; i < splines_.size(); i++) {
      state.position()[i] = splines_[i].getPositionAt(time);
      state.velocity()[i] = splines_[i].getVelocityAt(time);
      state.acceleration()[i] = splines_[i].getAccelerationAt(time);
    }

    return state;
  }

  void
  LocalPath::setDuration(double duration)
  {
    duration_ = duration;
    for(unsigned int i = 0; i < splines_.size(); i++)
      splines_[i].synchronize(duration_);
  }
}
