using Plots
ENV["GKSwstype"] = "100"
function createAspRatioHist(data_dir)
    aspect_ratio_0 = sort(reinterpret(Float64, read(data_dir * "aspect_ratio_0.bin")), rev=true)[1:1000]
    aspect_ratio_1 = sort(reinterpret(Float64, read(data_dir * "aspect_ratio_1.bin")), rev=true)[1:1000]
    aspect_ratio_2 = sort(reinterpret(Float64, read(data_dir * "aspect_ratio_2.bin")), rev=true)[1:1000]
    aspect_ratio_3 = sort(reinterpret(Float64, read(data_dir * "aspect_ratio_3.bin")), rev=true)[1:1000]

    bin_edges = 10:0.5:50
    hist = histogram(aspect_ratio_0,xlabel="aspect ratio", ylabel="freq.", alpha=0., bins=bin_edges, label="initial")
    histogram!(hist, aspect_ratio_1, xlabel="aspect ratio", ylabel="freq.", alpha=0.5, bins=bin_edges, label="1st smoothed")
    histogram!(hist, aspect_ratio_2, xlabel="aspect ratio", ylabel="freq.", alpha=0.5, bins=bin_edges, label="bisected")
    histogram!(hist, aspect_ratio_3, xlabel="aspect ratio", ylabel="freq.", alpha=0.5, bins=bin_edges, label="2nd smoothed")
    histogram(hist)
    savefig(data_dir * "fig/asp_hist.png")
end

createAspRatioHist("../tmp/result/")