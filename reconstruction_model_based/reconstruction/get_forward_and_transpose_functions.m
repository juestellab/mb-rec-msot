function [forward, transpose] = get_forward_and_transpose_functions(model, useSir)
%GET_FORWARD_AND_TRANSPOSE obtain forward and transpose functions as it is used in all reconstruction functions

if useSir
    forward = @model.forward;
    transpose = @model.transpose;
else
    forward = @model.forwardWithoutSir;
    transpose = @model.transposeWithoutSir;
end

end

