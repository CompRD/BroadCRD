#!/bin/sh

TestAssert NO_HEADER=True
TestAssert NO_HEADER=True TEST=eq USE_STRING=True
TestAssert NO_HEADER=True TEST=eq
TestAssert NO_HEADER=True TEST=ne USE_STRING=True
TestAssert NO_HEADER=True TEST=ne
TestAssert NO_HEADER=True TEST=lt
TestAssert NO_HEADER=True TEST=gt
TestAssert NO_HEADER=True TEST=le
TestAssert NO_HEADER=True TEST=ge
