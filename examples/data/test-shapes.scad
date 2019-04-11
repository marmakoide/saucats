module L_shape() {
  difference() {
    cube([20, 30, 10]);
    cube([20, 40, 40], center = true);
  }
}

module L_shape_container() {
  difference() {
    translate([-10, -10, 0])
    cube([40, 50, 20]);
    L_shape();
  }
}

L_shape_container();