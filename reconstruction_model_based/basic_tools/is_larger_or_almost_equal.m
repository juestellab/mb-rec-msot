function a_is_larger_than_b_or_almost_equal = is_larger_or_almost_equal(a, b, tolerance)

a_is_larger_than_b_or_almost_equal =  (a > b || abs(a-b) < tolerance);
end

