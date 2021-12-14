function a_is_smaller_than_b_or_almost_equal = is_smaller_or_almost_equal(a, b, tolerance)

a_is_smaller_than_b_or_almost_equal =  (a < b || abs(a-b) < tolerance);
end

