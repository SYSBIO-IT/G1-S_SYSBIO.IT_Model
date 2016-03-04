function Loss = Loss_Hill_genes(theta)

global x_axis y_axis

K = theta(1);
n = theta(2);

Hill_loss = x_axis.^n./(K^n + x_axis.^n);

Loss = sum((y_axis-Hill_loss).^2);