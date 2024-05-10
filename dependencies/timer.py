# ------------------------------------------------------------------------------
# Note: This file was copied from the factoritall repository by M. Ekerå.
#
# For further details, please see: https://github.com/ekera/factoritall
#
# Reference to other files and documentation below pertain to the above GitHub
# repository; not to the regevnum GitHub repository.

from time import time;
from math import floor;

# Supporting class to collect timing information.
class Timer:
  # A constant for state management indicating that the timer is stopped.
  STOPPED = 0;

  # A constant for state management indicating that the timer is running.
  RUNNING = 1;

  # Initializes the timer to a specific time delta, that defaults to zero.
  # The timer is left in the stopped state until manually started.
  def __init__(self, delta_t = 0):
    self.state = Timer.STOPPED;
    self.delta_t = delta_t;

  # Starts the timer if it is currently stopped. Returns the timer itself.
  def start(self):
    if self.state == Timer.STOPPED:
      self.state = Timer.RUNNING;
      self.t = time();

    return self;

  # Stops the timer if it is currently running. Returns the timer itself.
  def stop(self):
    if self.state == Timer.RUNNING:
      self.state = Timer.STOPPED;
      self.delta_t += time() - self.t;

    return self;

  # Stops the timer and re-initializes it to a specific time delta, that
  # defaults to zero. Returns the timer itself.
  def reset(self, delta_t = 0):
    self.state = Timer.STOPPED;
    self.delta_t = delta_t;

    return self;

  # Resets the timer and then starts it again. Returns the timer itself.
  def restart(self):
    self.reset();
    self.start();

    return self;

  # Peeks at the timer, returning the number of seconds elapsed. If the timer is
  # stopped, the time delta is returned. Otherwise, the sum of the time delta
  # and the current offset of the running timer is returned.
  def peek(self):
    tmp_delta_t = self.delta_t;
    if self.state == Timer.RUNNING:
      tmp_delta_t += time() - self.t;

    return tmp_delta_t;

  # Adds the time deltas of two stopped timers, returning a new timer with said
  # time delta. The timer is left in the stopped state until manually started.
  # If either timer is running, an exception is raised.
  def __add__(a, b):
    if (Timer.STOPPED != a.state) or (Timer.STOPPED != b.state):
      raise Exception("Error: Cannot add running timers.");

    return Timer(a.delta_t + b.delta_t);

  # Represents the timer as a string.
  def __repr__(self):
    # Get a temporary time delta.
    tmp_delta_t = self.peek();

    # Compute hours, minutes, seconds, milliseconds and microseconds.
    hours = floor(tmp_delta_t / 3600);
    mins = floor(tmp_delta_t / 60) % 60;
    secs = floor(tmp_delta_t) % 60;
    ms = int(floor((10 ** 3) * tmp_delta_t)) % (10 ** 3);
    us = int(floor((10 ** 6) * tmp_delta_t)) % (10 ** 3);

    # Format as a human-readable string.
    hr = "";
    if tmp_delta_t >= 3600:
      hr += str(hours) + " hour";
      if hours > 1:
        hr += "s ";
      else:
        hr += " ";

    if tmp_delta_t >= 60:
      hr += str(mins) + " min ";

    if tmp_delta_t >= 1:
      hr += str(secs) + " sec ";

    if tmp_delta_t >= 10 ** -3:
      hr += str(ms) + " ms ";

    hr += str(us) + " µs";

    return hr;
