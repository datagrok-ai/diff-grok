import {mrt, ros3prw, ros34prw} from '../../index';

export const methods = new Map([
  ['MRT', mrt],
  ['ROS3PRw', ros3prw],
  ['ROS34PRw', ros34prw],
]);

export const MAX_MAD = 0.1;
export const TIMEOUT_MS = 10000;
