
fprintf(format_number('F', -4, 'N'));       % Output: F = 42.00000
fprintf(format_number('F', -0.000123, 'N')); % Output: F = 123.00000u
fprintf(format_number('F', -150000000000000000000000000000000000000008, 'N')); % Output: F = 12.34568M

function formatted_string = format_number(variable_name, number, unit)
    
    if number == 0
        % No suffix for 0
        formatted_string = sprintf('%s = %10.5f %s\n', variable_name, number, unit);
    else
        % Extract the sign of the number
        sign = true;
        if number < 0
            sign = false;
            number = abs(number);
        end
    
        % Define small and big suffixes
        small_suffixes = {'', 'm', 'μ', 'n', 'p', 'f', 'a', 'z', 'y'};
        big_suffixes = {'', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};
    
        % Find the appropriate suffix and scale factor
        exp_value = floor(log10(number));
        idx_suffix = min(floor(abs(exp_value) / 3) + 1 + ((number < 1) * (1 - (mod(exp_value, 3) == 0))), 9);
        scale_factor = 10^( (1-2*(number < 1)) * 3 * (idx_suffix - 1) );
    
        % Format the number
        formatted_number = sprintf('%10.5f', (2*sign-1) * number / scale_factor);

        % Construct the formatted string
        formatted_string = sprintf('%s = %s %s%s\n', variable_name, formatted_number, (number < 1) * small_suffixes{idx_suffix} + (number >= 1) * big_suffixes{idx_suffix}, unit);
    end
end

