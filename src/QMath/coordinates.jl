function spherical_to_cartesian(r,θ,ϕ)
    x = r*sin(θ)*cos(ϕ)
    y = r*sin(θ)*sin(ϕ)
    z = r*cos(θ)

    [x, y, z]
end
