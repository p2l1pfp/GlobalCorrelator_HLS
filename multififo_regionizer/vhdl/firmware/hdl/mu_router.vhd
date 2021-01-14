library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity mu_router is
    generic(
            ETA_CENTER : integer 
    );
    port(
            ap_clk : IN STD_LOGIC;
            enabled : IN STD_LOGIC;
            newevent : IN STD_LOGIC;
            links_in      : IN glbparticles(NMUFIBERS-1 downto 0);
            fifo_in       : OUT particles(NPFREGIONS*NMUFIBERS-1 downto 0);
            fifo_in_write : OUT std_logic_vector(NPFREGIONS*NMUFIBERS-1 downto 0);
            fifo_in_roll  : OUT std_logic_vector(NPFREGIONS*NMUFIBERS-1 downto 0)
    );
end mu_router;

architecture Behavioral of mu_router is
    subtype wphi_t is signed(11 downto 0);
    subtype weta_t is signed(12 downto 0);
    constant ETA_SHIFT : weta_t := to_signed(-ETA_CENTER, weta_t'length);
begin

    link2fifo : process(ap_clk)
        variable local_phi : wphi_t;
        variable local_eta : weta_t;
    begin
        if rising_edge(ap_clk) then
            for ireg in 0 to NPFREGIONS-1 loop
                for ifib in 0 to NMUFIBERS-1 loop
                    -- convert to local coordinates
                    local_phi := links_in(ifib).phi - to_signed(ireg*PHI_SHIFT_INT, local_phi'length);
                    if local_phi < PHI_MPI then
                        local_phi := local_phi + PHI_2PI;
                    end if;
                    local_eta := links_in(ifib).eta + ETA_SHIFT;
                    -- prepare converted output
                    fifo_in(ireg*NMUFIBERS+ifib).pt   <= links_in(ifib).pt;
                    fifo_in(ireg*NMUFIBERS+ifib).eta  <= local_eta(9 downto 0);
                    fifo_in(ireg*NMUFIBERS+ifib).phi  <= local_phi(9 downto 0);
                    ---- CAREFUL with the rest: the good part is the high bits
                    fifo_in(ireg*NMUFIBERS+ifib).rest(27 downto 3) <= links_in(ifib).rest; 
                    fifo_in(ireg*NMUFIBERS+ifib).rest( 2 downto 0) <= (others => '0');
                    fifo_in_roll(ireg*NMUFIBERS+ifib) <= newevent;
                    if enabled = '1' and 
                           links_in(ifib).pt /= 0 and
                           local_phi <= PHI_HALFWIDTH_POS and
                           local_phi >= PHI_HALFWIDTH_NEG and
                           local_eta <= ETA_HALFWIDTH_POS and
                           local_eta >= ETA_HALFWIDTH_NEG then
                        fifo_in_write(ireg*NMUFIBERS+ifib) <= '1';
                    else
                        fifo_in_write(ireg*NMUFIBERS+ifib) <= '0';
                    end if;
                end loop; 
            end loop; 
        end if;
   end process link2fifo;


end Behavioral;
