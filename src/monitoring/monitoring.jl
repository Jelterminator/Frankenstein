module Monitoring

using ..Core
export NullMonitor

struct NullMonitor <: AbstractPerformanceMonitor end

end
