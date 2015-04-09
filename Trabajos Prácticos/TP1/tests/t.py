#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

def main():
  f = open('test1.in', "rU")
  lists = []
  n = 0
  nl = f.readline()
  for line in f:
    for number in line.split():
      n += 1
      if n == 4:
        if not(float(number) in lists):
          lists.append(float(number))
        n = 0

  f.close()
  print 'Sanguis:',
  print lists

  print 'Expected'
  d = open('test1.expected','rU')
  n = 0
  final = []
  for line in d:
    l = line.split()
    if float(l[2]) in lists:
      final.append([float(l[0]),float(l[1]),float(l[2])])
      print 'X:', l[0], 'Y:', l[1], 'Temp:', l[2]
  d.close()

  print 'Out'
  d = open('test1.out','rU')
  n = 0
  for line in d:
    l = line.split()
    if float(l[2]) in lists:
      print 'X:', l[0], 'Y:', l[1], 'Temp:', l[2]
  d.close()

if __name__ == '__main__':
  main()
