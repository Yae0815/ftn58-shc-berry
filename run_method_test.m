% Test batch script to verify Bastin vs Weighted methods
clear;
clc;

fprintf('=== Testing Bastin vs Weighted Methods ===\n');
fprintf('This will test if the two methods give different results\n');
fprintf('at finite temperature (T=300K)\n\n');

% Run the test
test_bastin_vs_weighted();

fprintf('\n=== Test Complete ===\n');