
begin
    count = 0
    s = Set{Int}()
    N = 4
    for t in Iterators.product([0:N for _ in 1:10]...)
        if sum(t) <= N
            count += 1
            u = (t[1] | (t[2] << 3) | (t[3] << 6) | (t[4] << 9) | (t[5] << 12) | 
                (t[6] << 15) | (t[7] << 18) | (t[8] << 21) | (t[9] << 24) | (t[10] << 27))
            push!(s, u)
        end
    end
    count, length(s)
end

N_PIECES = 4
binomial(10+(N_PIECES-1), N_PIECES)

binomial(10, 3)

l = collect(s)

import Random
begin
    Random.seed!(0)
    p = 2029
    best_n = 0
    best_k = 0
    for k in 1:(2^16)
        # k = rand(Int)
        n_keys = length(Set([(e ⊻ k) % p for e in l]))
        if n_keys > best_n
            best_n = n_keys
            best_k = k
        end
    end
    best_n, best_k
end
hashes = sort([(e ⊻ best_k) % p for e in l])
for i in 1:length(hashes)-1
    if hashes[i] == hashes[i+1]
        n_dups = 0
        i_end = i
        for j in i:length(hashes)-1
            if hashes[i] == hashes[j]
                n_dups += 1
            else
                i_end = j
                break
            end
        end
        println(i, " ", hashes[i], " ", n_dups, " -> ", hashes[i_end])
    end
end

# 16 workers unpinned: 2090s
# 16 workers pinned: 2073s
# 32 workers pinned: 2195s